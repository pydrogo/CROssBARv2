#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - CROssBAR prototype
"""

import pandas as pd
import numpy as np
from tqdm import tqdm  # Progress bar

import os
import yaml
import importlib as imp

import biocypher
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class BiocypherAdapter:
    """
    The connection can be defined in three ways:
        * Providing a ready ``neo4j.Driver`` instance
        * By URI and authentication data
        * By a YML config file

    Args:
        driver (neo4j.Driver): A ``neo4j.Driver`` instance, created by,
            for example, ``neo4j.GraphDatabase.driver``.
        db_name (str): Name of the database (Neo4j graph) to use.
        db_uri (str): Protocol, host and port to access the Neo4j server.
        db_auth (tuple): Neo4j server authentication data: tuple of user
            name and password.
        config_file (str): Path to a YML config file which provides the URI,
            user name and password.
        network (pypath.core.network.Network): A network database object.
        wipe (bool): Wipe the database after connection, ensuring the data
            is loaded into an empty database.
    """

    def __init__(
        self,
        driver=None,
        db_name="neo4j",
        db_uri="bolt://localhost:7687",
        db_user="neo4j",
        db_passwd="password_here",
        user_schema_config_path= r"config/schema_config.yaml",
        wipe=True,
        offline=False,
        delimeter = ";",
        array_delimiter = "|",
        quote_char = "'",
        skip_bad_relationships = True,
        skip_duplicate_nodes = True,
    ):

        self.bcy = biocypher.Driver(
            driver=driver,
            db_name=db_name,
            db_uri=db_uri,
            db_user=db_user,
            db_passwd=db_passwd,
            user_schema_config_path=user_schema_config_path,
            wipe=wipe,
            offline=offline,
            delimiter=delimeter,
            array_delimiter=array_delimiter,
            quote_char=quote_char,
            skip_bad_relationships=skip_bad_relationships,
            skip_duplicate_nodes=skip_duplicate_nodes,
        )

        self.nodes = None
        self.edges = None

    def build_python_object(self, data=None, early_stopping=None):
        """
        Loads uniprot data and preprocess it
        """
        protein_nodes, gene_nodes, organism_nodes, gene_to_protein_edges, protein_to_organism_edges = generate_nodes_and_edges(data, early_stopping=early_stopping)
        
        self.nodes = protein_nodes + gene_nodes + organism_nodes
        self.edges = gene_to_protein_edges + protein_to_organism_edges

    def translate_python_object_to_neo4j(self, nodes=None, edges=None):
        """
        Loads a pypath network into the biocypher (Neo4j) backend.

        Args:
            - network (pypath.core.network.Network): A network database
              object. If `None`, the value of :py:attr:`network` will be
              used.
        """
        if nodes is None:
            nodes = self.nodes
        if edges is None:
            edges = self.edges


    def write_nodes(self, nodes=None, db_name="neo4j"):
        """
        Writes biocypher nodes to CSV files for admin import.

        Args:
            nodes (list): A list of nodes.
            db_name (str): Name of the database (Neo4j graph) to use.
        """

        if nodes is None:
            nodes = self.nodes

        id_type_tuples = list(process_nodes(nodes))

        # write nodes
        self.bcy.write_nodes(id_type_tuples, db_name=db_name)

    def write_edges(self, edges=None, db_name="neo4j"):
        """
        Writes biocypher edges to CSV files for admin import.

        Args:
            edges (list): A list of edges.
            db_name (str): Name of the database (Neo4j graph) to use.
        """
        if edges is None:
            edges = self.edges

        src_tar_type_tuples = list(process_edges(edges))

        self.bcy.write_edges(src_tar_type_tuples, db_name=db_name)

    def write_to_csv_for_admin_import(
        self, nodes=None, edges=None, db_name="neo4j"
    ):
        """
        Loads a pypath network into the biocypher (Neo4j) backend using
        the fast Admin Import function, which requires text files that
        need to be properly formatted since it turns off safety measures
        at import.

        Args:
            nodes (list): A list of nodes.
            edges (list): A list of edges.
            db_name (str): Name of the database (Neo4j graph) to use.
        """
        if nodes is None:
            nodes = self.nodes

        if edges is None:
            edges = self.edges

        self.write_nodes(nodes, db_name)
        self.write_edges(edges, db_name)
        self.bcy.write_import_call()

    def load(self, obj):
        """
        Loads any compatible object into the biocypher (Neo4j) database.

        Args:
            obj: An object from this module compatible with the current
                adapter. Currently the following database objects are
                supported:
                    * :py:class:`pypath.core.network.Network`
        """

        if hasattr(obj, "nodes") and hasattr(obj, "interactions"):

            self.translate_python_object_to_neo4j(network=obj)



def process_nodes(nodes):
    for node_dict in nodes:
        if "accession" in node_dict.keys():
            _id = str(node_dict["accession"])
            _type = "Protein"
            _props = {}
            for k, v in node_dict.items():
                if k == "accession" or str(v) == "nan":
                    continue
                else:
                    if "|" not in str(v):
                        _props[k] = v
                    else:
                        _props[k] = str(v).split("|")
            
        elif "entrez_id" in node_dict.keys():
            _id = str(node_dict["entrez_id"])
            _type = "Gene"
            _props = {}
            for k, v in node_dict.items():
                if k == "entrez_id" or str(v) == "nan":
                    continue
                else:
                    if "|" not in str(v):
                        _props[k] = v
                    else:
                        _props[k] = str(v).split("|")
        else:
            _id = str(node_dict["tax_id"])
            _type = "Organism"
            _props = {"organism":node_dict["organism"]}
                        
        yield (_id, _type, _props)

def process_edges(edges):
    _props = {}
    for edge_dict in edges:
        if "entrez_id" in edge_dict.keys():
            _source = str(edge_dict["entrez_id"])
            _target = str(edge_dict["accession"])
            _type = "Encodes"
        else:
            _source = str(edge_dict["accession"])
            _target = str(edge_dict["tax_id"])
            _type = "Belongs_To"
        
        yield (_source, _target, _type, _props)


def generate_nodes_and_edges(uniprot_df, early_stopping=None):
    node_values = ["accession", "secondary_accessions", "length", "mass", "tax_id", "organism", "protein_names", "chromosome_ids", 
                   "chromosome_info", "entrez_id", "ec_numbers", "kegg", "ensembl_transcript", "ensembl_gene_ids", "genes", 
                   "virus_hosts_tax_ids"] # list of attributes that will be added graph database

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
        entrez_id, ec_numbers, kegg, ensembl_transcript, ensemble_gene_ids, genes, virus_hosts_tax_ids = clean_uniprot_data_for_nodes(row)

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

        # replace sensitive elements for admin-import csv from protein nodes
        for k, v in protein_add_dict.items():
            if isinstance(v, list):
                if len(v) == 1:
                    protein_add_dict[k] = str(v[0]).replace(";",":").replace("|",",").replace("'","")
                else:
                    v_update = [str(e).replace("|",",").replace("'","") for e in v]
                    protein_add_dict[k] = "|".join(v_update).replace(";",":")

            elif isinstance(v, str):
                protein_add_dict[k] = str(v).replace(";",":").replace("|",",").replace("'","")


        # replace sensitive elements for admin-import csv from gene nodes
        for k, v in gene_add_dict.items():
            if isinstance(v, list):
                if len(v) == 1:
                    gene_add_dict[k] = str(v[0]).replace(";",":").replace("|",",").replace("'","")
                else:
                    v_update = [str(e).replace("|",",").replace("'","") for e in v]
                    gene_add_dict[k] = "|".join(v_update).replace(";",":")

            elif isinstance(v, str):
                gene_add_dict[k] = str(v).replace(";",":").replace("|",",").replace("'","")

        # replace sensitive elements for admin-import csv from organism nodes
        for k, v in organism_add_dict.items():
            organism_add_dict[k] = str(v).replace("'","")
            
        # create protein-organism edges
        protein_to_organism_interaction_dict = {"accession":accession, "tax_id":tax_id}
        protein_to_organism_interactions.append(protein_to_organism_interaction_dict)


        # check duplicates in protein node, if not add it to list
        if protein_add_dict not in protein_check_duplicates_admin_import:
            protein_check_duplicates_admin_import.append(protein_add_dict)

        # check genes and entrez_id is whether NaN and check gene node whether duplicate, if not create gene node
        if str(genes) != "nan" and str(entrez_id) != "nan" and gene_add_dict not in gene_check_duplicates_admin_import and entrez_id not in entrez_id_check_duplicates:
            entrez_id_check_duplicates.append(entrez_id)
            gene_check_duplicates_admin_import.append(gene_add_dict)

        # create gene-protein edges
        if str(genes) != "nan" and str(entrez_id) != "nan":
            gene_to_protein_interaction_dict = {"entrez_id":entrez_id, "accession":accession}
            gene_to_protein_interactions.append(gene_to_protein_interaction_dict)

        # create organism nodes if not duplicate
        if organism_add_dict not in organism_check_duplicates_admin_import:
            organism_check_duplicates_admin_import.append(organism_add_dict)
    
    return protein_check_duplicates_admin_import, gene_check_duplicates_admin_import, organism_check_duplicates_admin_import,\
           gene_to_protein_interactions, protein_to_organism_interactions

def clean_uniprot_data_for_nodes(row):
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
    # CURRENT FORMAT:INT
    # returns length
    length = int(row["length"]) 
    
    # mass
    # CURRENT FORMAT:INT
    # returns mass
    mass = int(row["mass"]) 
    
    # tax_id
    # CURRENT FORMAT:INT
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
        
    
    # if there are multiple ensembl gene ids reduce them with Opentargets
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
