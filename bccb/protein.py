import os
import sys
import collections
import pandas as pd
from time import time
from pathlib import Path

from pypath.inputs import uniprot
from pypath.utils import mapping
from pypath.share import curl
import biocypher

from adapter import BiocypherAdapter

from tqdm import tqdm # progress bar


class Uniprot_data:
    def __init__(self, organism='*', rev = True):        
        self.organism = organism
        self.rev = rev
        self.uniprot_df = None
    
    def xref_process(self, attribute_key, protein):
    # delete trailing ';' from fields containing only one id
        xrefs = self.data[attribute_key].get(protein)
        if xrefs:
            xrefs_splitted = xrefs.split(';')
            xrefs_splitted = list(filter(None, xrefs_splitted))
            if len(xrefs_splitted) > 1:
                return xrefs
            else:
                return xrefs_splitted[0]


    
    def ensembl_process(self, enst_list): 
    # take ensembl transcript ids, return ensembl gene ids by using pypath mapping tool
        ensg_ids = set()
        for enst_id in enst_list.split(';'):
            only_transcript_id = enst_id.split(' [')[0] # strip from uniprot alternative transcript id
            only_transcript_id = only_transcript_id.split('.')[0]
            ensg_id= list(mapping.map_name(only_transcript_id, 'enst_biomart', 'ensg_biomart'))
            ensg_id= (
                ensg_id[0]
                    if ensg_id else
                None
            )
            if ensg_id:
                ensg_ids.add(ensg_id)
        
        ensg_list = (
            ';'.join(ensg_ids)
                if ensg_ids else
            None
        )

        return ensg_list


    def uniprot_data_download(self):
        
        t0 = time()

        self.attributes = [ # with query field keys in uniprot api
            # primary attributes
            'length', 'mass', 'organism','organism-id', 'protein names', 'proteome', 'genes', 'ec', 'database(Ensembl)',
            # xref attributes
            'database(GeneID)', 'virus hosts','database(OpenTargets)', 'database(KEGG)', 'database(DisGeNET)',
            'database(EnsemblBacteria)', 'database(EnsemblFungi)', 'database(EnsemblMetazoa)',
            'database(EnsemblPlants)', 'database(EnsemblProtists)',
        ]

        # download all swissprot ids
        with curl.cache_off():
            self.uniprot_ids=list(uniprot._all_uniprots(self.organism, self.rev))
        
        # download attribute dicts
        self.data = {}
        for query_key in tqdm(self.attributes):
            try: # try downloading two times in case of connection issues etc
                with curl.cache_off():
                    self.data[query_key] = uniprot.uniprot_data(query_key, self.organism, self.rev)
            except Exception as e:
                with curl.cache_off():
                    self.data[query_key] = uniprot.uniprot_data(query_key, self.organism, self.rev)
            print(f' {query_key} field is downloaded')
        
        secondary_ids = uniprot.get_uniprot_sec(None)
        self.data['secondary_ids'] = collections.defaultdict(list)
        for id in secondary_ids:
            self.data['secondary_ids'][id[1]].append(id[0])
        for k,v in self.data['secondary_ids'].items():
            self.data['secondary_ids'][k] = ';'.join(v)

        t1 = time()
        print(f' Downloaded data from UniProtKB in {round((t1-t0) / 60, 2)} mins, now writing csv files')

    
    def build_dataframe(self):
        print("Building dataframe")
        protein_dict_list = []
        for protein in tqdm(self.uniprot_ids):
            protein_dict = {
                    # define primary attributes which needs specific preprocessing
                    'accession': protein,
                    'secondary_accessions': self.data['secondary_ids'].get(protein),
                    'length': int(self.data['length'].get(protein)),
                    'mass': (
                            int(self.data['mass'][protein].replace(',',''))
                                if protein in self.data['mass'] else
                            None
                        ),
                    'tax_id': int(self.data['organism-id'].get(protein)),
                    'organism': self.data['organism'].get(protein),
                    'protein_names': self.data['protein names'].get(protein),
                    'chromosome': (
                            ';'.join(self.data['proteome'].get(protein).split(','))
                                if self.data['proteome'].get(protein) else
                            None
                        ),
                    'genes': (
                            ';'.join(self.data['genes'][protein].split(' '))
                                if protein in self.data['genes'] else
                            None
                        ),

                    'ec_numbers': self.data['ec'].get(protein),
                    'ensembl_transcript': self.xref_process('database(Ensembl)', protein),
                    'ensembl_gene': (
                            self.ensembl_process(self.xref_process('database(Ensembl)', protein))
                                if self.xref_process('database(Ensembl)', protein) else
                            None
                        ),
                }
            
            protein_dict_list.append(protein_dict)
        
        self.uniprot_df = pd.DataFrame.from_dict(protein_dict_list, orient= 'columns') # create uniprot dataframe
        
        self.uniprot_df.replace("nan", np.nan, inplace=True)
        self.uniprot_df.fillna(value=np.nan, inplace=True) # replace None with NaN
        print("Dataframe is builded")
        
    def generate_nodes_and_edges(self, uniprot_df=None, early_stopping=None):
        
        if uniprot_df is None:
            uniprot_df = self.uniprot_df
        
        node_values = ["accession", "secondary_accessions", "length", "mass", "tax_id", "organism", "protein_names", "chromosome_ids", 
                       "chromosome_info", "entrez_id", "ec_numbers", "kegg", "ensembl_transcript", "ensembl_gene_ids", "genes", 
                       "virus_hosts_tax_ids"] # list of values that will be added graph database

        keys = list(range(0, len(node_values)))
        order = dict()
        for k, v in zip(keys, node_values):
            order[k] = v

        # attributes of nodes (protein, gene, organism)
        protein_attributes_list = ["accession", "secondary_accessions", "length", "mass", "protein_names", "chromosome_ids", 
                                   "chromosome_info", "tax_id", "ec_numbers", "virus_hosts_tax_ids"] 
        organism_attributes_list = ["tax_id", "organism"]
        gene_attributes_list = ["entrez_id", "genes", "kegg", "ensembl_transcript", "ensembl_gene_ids", "tax_id"]


        # lists for checking duplicates and creating nodes
        protein_check_duplicates_admin_import = []
        gene_check_duplicates_admin_import = []
        organism_check_duplicates_admin_import = []
        entrez_id_check_duplicates = []

        # lists for creating edges
        gene_to_protein_interactions = []
        protein_to_organism_interactions = []

        for index, row in tqdm(uniprot_df.iterrows()):

            if early_stopping is not None and index == early_stopping:
                break

            accession, secondary_accessions, length, mass, tax_id, organism, protein_names, chromosome_ids, chromosome_info,\
            entrez_id, ec_numbers, kegg, ensembl_transcript, ensemble_gene_ids, genes, virus_hosts_tax_ids = self.clean_uniprot_data_for_nodes(row)

            all_values = [accession, secondary_accessions, length, mass, tax_id, organism, protein_names, chromosome_ids, 
                          chromosome_info, entrez_id, ec_numbers, kegg, ensembl_transcript, ensemble_gene_ids, genes, 
                          virus_hosts_tax_ids]

            # check if there is NaN, if not add them belonging nodes (i.e., protein, gene, organism)
            protein_add_dict = {}
            gene_add_dict = {}
            organism_add_dict = {}
            for idx, v in enumerate(all_values):        
                if order[idx] in protein_attributes_list:                
                    protein_add_dict[order[idx]] = v
                if order[idx] in gene_attributes_list:
                    gene_add_dict[order[idx]] = v
                if order[idx] in organism_attributes_list:
                    organism_add_dict[order[idx]] = v

            # replace sensitive elements for admin-import csv from protein node
            for k, v in protein_add_dict.items():
                if isinstance(v, list):
                    if len(v) == 1:
                        protein_add_dict[k] = str(v[0]).replace(";",":").replace("|",",").replace("'","")
                    else:
                        v_update = [str(e).replace("|",",").replace("'","") for e in v]
                        protein_add_dict[k] = "|".join(v_update).replace(";",":")

                elif isinstance(v, str):
                    protein_add_dict[k] = str(v).replace(";",":").replace("|",",").replace("'","")


            # replace sensitive elements for admin-import csv from gene node
            for k, v in gene_add_dict.items():
                if isinstance(v, list):
                    if len(v) == 1:
                        gene_add_dict[k] = str(v[0]).replace(";",":").replace("|",",").replace("'","")
                    else:
                        v_update = [str(e).replace("|",",").replace("'","") for e in v]
                        gene_add_dict[k] = "|".join(v_update).replace(";",":")

                elif isinstance(v, str):
                    gene_add_dict[k] = str(v).replace(";",":").replace("|",",").replace("'","")

            # replace sensitive elements for admin-import csv from organism node
            for k, v in organism_add_dict.items():
                organism_add_dict[k] = str(v).replace("'","")

            # create protein-organism edge
            protein_to_organism_interaction_dict = {"accession":accession, "tax_id":tax_id}
            protein_to_organism_interactions.append(protein_to_organism_interaction_dict)


            # check duplicates in protein node, if not add it to list
            if protein_add_dict not in protein_check_duplicates_admin_import:
                protein_check_duplicates_admin_import.append(protein_add_dict)

            # check genes and entrez_id is whether NaN and check gene node whether duplicate, if not create gene node
            if str(genes) != "nan" and str(entrez_id) != "nan" and gene_add_dict not in gene_check_duplicates_admin_import and entrez_id not in entrez_id_check_duplicates:
                entrez_id_check_duplicates.append(entrez_id)
                gene_check_duplicates_admin_import.append(gene_add_dict)

            # create gene-protein edge
            if str(genes) != "nan" and str(entrez_id) != "nan":
                gene_to_protein_interaction_dict = {"entrez_id":entrez_id, "accession":accession}
                gene_to_protein_interactions.append(gene_to_protein_interaction_dict)

            # create organism node if not duplicate
            if organism_add_dict not in organism_check_duplicates_admin_import:
                organism_check_duplicates_admin_import.append(organism_add_dict)

        return protein_check_duplicates_admin_import, gene_check_duplicates_admin_import, organism_check_duplicates_admin_import,\
               gene_to_protein_interactions, protein_to_organism_interactions

    def clean_uniprot_data_for_nodes(self, row):
        # accession
        # CURRENT FORMAT:STR
        # returns accession
        accession = str(row["accession"])


        # secondary_accessions
        # CURRENT FORMAT:LIST
        # returns secondary_accessions
        if str(row["secondary_accessions"]) == "nan":
            secondary_accessions = np.nan
        elif ";" in str(row["secondary_accessions"]):
            secondary_accessions = str(row["secondary_accessions"]).strip().split(";")
        else:
            secondary_accessions = [str(row["secondary_accessions"]).strip()]

        # length
        # CURRENT FORMAT:INT (MAY CHANGE)
        # returns length
        length = int(row["length"]) 

        # mass
        # CURRENT FORMAT:INT (MAY CHANGE)
        # returns mass
        mass = int(row["mass"]) 

        # tax_id
        # CURRENT FORMAT:INT (MAY CHANGE)
        # returns tax_id
        tax_id = int(row["tax_id"]) 

        # organism 
        # CURRENT FORMAT:STR
        # returns organism
        if str(row["organism"]) == "nan":
            organism = np.nan
        else:
            organism = str(row["organism"])


        # protein_names
        # CURRENT FORMAT:LIST
        # returns protein_names
        if "[Cleaved" in str(row["protein_names"]):
            # discarding part after the "[Cleaved"
            clip_index = str(row["protein_names"]).index("[Cleaved")
            protein_names = [str(row["protein_names"])[:clip_index].replace("(Fragment)","").strip()]

            # handling multiple protein names
            if "(EC" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []

                for name in splitted:                    
                    if not name.strip().startswith("EC"):
                        if not name.strip().startswith("Fragm"):
                            if name.strip().endswith(")"):
                                protein_names.append(name.replace(")","").strip())
                            else:
                                protein_names.append(name.strip())

            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []
                for name in splitted:
                    if not name.strip().startswith("Fragm"):
                        if name.strip().endswith(")"):
                            protein_names.append(name.replace(")","").strip())
                        else:
                            protein_names.append(name.strip())

        elif "[Includes" in str(row["protein_names"]):
            # discarding part after the "[Includes"
            clip_index = str(row["protein_names"]).index("[Includes")
            protein_names = [str(row["protein_names"])[:clip_index].replace("(Fragment)","").strip()]
            # handling multiple protein names
            if "(EC" in protein_names[0]:

                splitted = protein_names[0].split(" (")
                protein_names = []

                for name in splitted:                    
                    if not name.strip().startswith("EC"):
                        if not name.strip().startswith("Fragm"):
                            if name.strip().endswith(")"):                            
                                protein_names.append(name.replace(")","").strip())
                            else:
                                protein_names.append(name.strip())

            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []
                for name in splitted:
                    if not name.strip().startswith("Fragm"):
                        if name.strip().endswith(")"):                        
                            protein_names.append(name.replace(")","").strip())
                        else:
                            protein_names.append(name.strip())



        # handling multiple protein names
        elif "(EC" in str(row["protein_names"]).replace("(Fragment)",""):
            splitted = str(row["protein_names"]).split(" (")
            protein_names = []

            for name in splitted:                
                if not name.strip().startswith("EC"):
                    if not name.strip().startswith("Fragm"):
                        if name.strip().endswith(")"):
                            protein_names.append(name.replace(")","").strip())
                        else:
                            protein_names.append(name.strip())

        elif " (" in str(row["protein_names"]).replace("(Fragment)",""):
            splitted = str(row["protein_names"]).split(" (")
            protein_names = []
            for name in splitted:
                if not name.strip().startswith("Fragm"):
                    if name.strip().endswith(")"):                    
                        protein_names.append(name.replace(")","").strip())
                    else:
                        protein_names.append(name.strip())


        else:
            protein_names = [str(row["protein_names"]).replace("(Fragment)","").strip()]


        # chromosome
        # CURRENT FORMAT:LIST
        # returns chromosome_ids, chromosome_info
        if str(row["chromosome"]) == "nan":
            chromosome_ids = np.nan
            chromosome_info = np.nan

        elif ";" in str(row["chromosome"]):
            ch_list = str(row["chromosome"]).split(";")
            chromosome_ids = []
            chromosome_info = []
            for pair in ch_list:
                if len(pair.split(":")) == 2:
                    chromosome_ids.append(pair.split(":")[0].strip())
                    chromosome_info.append(pair.split(":")[1].strip())

            if chromosome_info.count("Chromosome") > 1:
                chromosome_info = list(set(chromosome_info))

        else:
            chromosome_ids = [str(row["chromosome"]).split(":")[0].strip()]
            chromosome_info = [str(row["chromosome"]).split(":")[1].strip()]


        # genes
        # CURRENT FORMAT:LIST
        # returns genes
        if str(row["genes"]) == "nan":
            genes = np.nan
        elif ";" in str(row["genes"]):
            genes = str(row["genes"]).split(";")
        else:
            genes = [str(row["genes"])]


        # ec_numbers
        # CURRENT FORMAT:LIST
        # returns ec_numbers
        if str(row["ec_numbers"]) == "nan":
            ec_numbers = np.nan        
        elif ";" in str(row["ec_numbers"]):
            ec_numbers = [e.strip() for e in str(row["ec_numbers"]).split(";")]
        else:
            ec_numbers = [str(row["ec_numbers"]).strip()]


        # kegg
        # CURRENT FORMAT:LIST
        # returns kegg
        if str(row["database(KEGG)"]) == "nan":
            kegg = np.nan
        elif ";" in str(row["database(KEGG)"]):
            splitted = str(row["database(KEGG)"])[:-1].split(";")
            kegg = []
            for pair in splitted:
                if str(genes) == "nan":
                    kegg.append(pair.split(":")[1].strip())
                elif pair.split(":")[1].strip() in genes:
                    kegg.append(pair.split(":")[1].strip())
                else:
                    kegg.append(pair.split(":")[1].strip())
        else:
            kegg = [str(row["database(KEGG)"]).split(":")[1].strip()]


        # ensembl_transcript
        # CURRENT FORMAT:LIST
        # returns ensembl_transcript
        if str(row["ensembl_transcript"]) == "nan":
            ensembl_transcript = np.nan
        elif ";" in str(row["ensembl_transcript"]):
            splitted = str(row["ensembl_transcript"])[:-1].split(";")
            ensembl_transcript = [tr.strip() for tr in splitted]
        else:
            ensembl_transcript = [str(row["ensembl_transcript"]).strip()]


        # ensembl_gene
        # CURRENT FORMAT:LIST
        # returns ensembl_gene
        if str(row["ensembl_gene"]) == "nan":
            ensembl_gene = np.nan
        elif ";" in str(row["ensembl_gene"]):
            splitted = str(row["ensembl_gene"]).split(";")
            ensembl_gene = [g.strip() for g in splitted]
        else:
            ensembl_gene = [str(row["ensembl_gene"]).strip()]


        # entrez_id
        # CURRENT FORMAT:STR
        # returns entrez_id
        if str(row["database(GeneID)"]) == "nan":
            entrez_id = np.nan
        elif ";" in str(row["database(GeneID)"]):
            splitted = str(row["database(GeneID)"])[:-1].split(";")
            entrez_id = [ent.strip() for ent in splitted]
            if str(row["database(DisGeNET)"]) != "nan" and ";" not in str(row["database(DisGeNET)"]) and str(row["database(DisGeNET)"]) in entrez_id:
                entrez_id = str(row["database(DisGeNET)"]) # takes the DisGeNET's entrez_id as id 
            else:            
                entrez_id = str(entrez_id[0]) # takes first element as entrez_id
        else:
            entrez_id = str(row["database(GeneID)"]).strip()

        # EnsemblBacteria
        # CURRENT FORMAT:LIST
        # returns EnsemblBacteria
        if str(row["database(EnsemblBacteria)"]) == "nan":
            EnsemblBacteria = np.nan
        elif ";" in str(row["database(EnsemblBacteria)"]):
            splitted = str(row["database(EnsemblBacteria)"])[:-1].split(";")
            EnsemblBacteria = [ens.strip() for ens in splitted]
        else:
            EnsemblBacteria = [str(row["database(EnsemblBacteria)"]).strip()]

        # EnsemblFungi
        # CURRENT FORMAT: LIST
        # returns EnsemblFungi
        if str(row["database(EnsemblFungi)"]) == "nan":
            EnsemblFungi = np.nan
        elif ";" in str(row["database(EnsemblFungi)"]):
            splitted = str(row["database(EnsemblFungi)"])[:-1].split(";")
            EnsemblFungi = [ens.strip() for ens in splitted]
        else:
            EnsemblFungi = [str(row["database(EnsemblFungi)"]).strip()]

        # EnsemblMetazoa
        # CURRENT FORMAT: LIST
        # returns EnsemblMetazoa
        if str(row["database(EnsemblMetazoa)"]) == "nan":
            EnsemblMetazoa = np.nan
        elif ";" in str(row["database(EnsemblMetazoa)"]):
            splitted = str(row["database(EnsemblMetazoa)"])[:-1].split(";")
            EnsemblMetazoa = [ens.strip() for ens in splitted]
        else:
            EnsemblMetazoa = [str(row["database(EnsemblMetazoa)"]).strip()]

        # EnsemblPlants
        # CURRENT FORMAT: LIST
        # returns EnsemblPlants
        if str(row["database(EnsemblPlants)"]) == "nan":
            EnsemblPlants = np.nan
        elif ";" in str(row["database(EnsemblPlants)"]):
            splitted = str(row["database(EnsemblPlants)"])[:-1].split(";")
            EnsemblPlants = [ens.strip() for ens in splitted]
        else:
            EnsemblPlants = [str(row["database(EnsemblPlants)"]).strip()]


        # EnsemblProtists
        # CURRENT FORMAT: LIST
        # returns EnsemblProtists
        if str(row["database(EnsemblProtists)"]) == "nan":
            EnsemblProtists = np.nan
        elif ";" in str(row["database(EnsemblProtists)"]):
            splitted = str(row["database(EnsemblProtists)"])[:-1].split(";")
            EnsemblProtists = [ens.strip() for ens in splitted]
        else:
            EnsemblProtists = [str(row["database(EnsemblProtists)"]).strip()]


        # virus_hosts
        # CURRENT FORMAT:LIST
        # returns virus_hosts_tax_ids
        if str(row["virus hosts"]) == "nan":
            virus_hosts_tax_ids = np.nan
        elif ";" in str(row["virus hosts"]):
            splitted = str(row["virus hosts"]).split(";")
            virus_hosts_tax_ids = []
            for v in splitted:
                virus_hosts_tax_ids.append(v[v.index("[")+1:v.index("]")].split(":")[1].strip())
        else:
            virus_hosts = str(row["virus hosts"])
            virus_hosts_tax_ids = [virus_hosts[virus_hosts.index("[")+1:virus_hosts.index("]")].split(":")[1].strip()]


        # if there are multiple ensemble_gene ids reduce them with Opentargets
        if str(ensembl_gene) != "nan" and len(ensembl_gene) > 1 and str(row["database(OpenTargets)"]) != "nan" and ";" not in str(row["database(OpenTargets)"]) and str(row["database(OpenTargets)"]) in ensembl_gene:
            ensembl_gene = [str(row["database(OpenTargets)"])]


        # concat all ensemble gene ids
        ensemble_gene_ids = []
        concat_list = [ensembl_gene, EnsemblBacteria, EnsemblFungi, EnsemblMetazoa, EnsemblPlants, EnsemblProtists]
        for element in concat_list:
            if str(element) != "nan":
                ensemble_gene_ids.extend(element)

        if len(ensemble_gene_ids) == 0: # if nothing in there make it NaN
            ensemble_gene_ids = np.nan


        """
        print("accession:", accession)
        print("secondary_accessions:", secondary_accessions)
        print("length:", length)
        print("mass:", mass)
        print("tax_id:", tax_id)
        print("organism:", organism)
        print("before protein_names:", str(row["protein_names"]))
        print("protein_names:", protein_names)
        print("chromosome_ids:", chromosome_ids)
        print("chromosome_info:", chromosome_info)
        print("genes:", genes)
        print("ec_numbers:", ec_numbers)
        print("kegg:", kegg)
        print("ensembl_transcript:", ensembl_transcript)
        print("all ensemble ids:", ensemble_gene_ids)
        print("entrez_id:", entrez_id)
        print("before virus_hosts:", str(row["virus hosts"]))
        print("virus_hosts_tax_ids:", virus_hosts_tax_ids)
        """

        return accession, secondary_accessions, length, mass, tax_id, organism, protein_names, chromosome_ids, chromosome_info,\
               entrez_id, ec_numbers, kegg, ensembl_transcript, ensemble_gene_ids, genes, virus_hosts_tax_ids
    
    def build_nodes_and_edges(self, early_stopping=None):
        # generate nodes and edges
        protein_nodes, gene_nodes, organism_nodes, gene_to_protein_edges, protein_to_organism_edges = self.generate_nodes_and_edges(data, early_stopping=early_stopping)
        
        self.nodes = protein_nodes + gene_nodes + organism_nodes
        self.edges = gene_to_protein_edges + protein_to_organism_edges
        
    def call_biocypher_adapter(self):
        """
        if "early_stopping" is specified with an integer value, it stops preprocessing on the specified row.
        """

        print("Calling Biocypher adapter")
        adapt = BiocypherAdapter(offline=True, db_name="neo4j", wipe=True, quote_char="'")
        adapt.write_to_csv_for_admin_import(nodes=self.nodes, edges=self.edges)
