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
    def call_biocypher_adapter(self, early_stopping=None):
        """
        if "early_stopping" is specified with an integer value, it stops preprocessing on the specified row.
        """

        print("Calling Biocypher adapter")
        adapt = BiocypherAdapter(offline=True, db_name="neo4j", wipe=True, quote_char="'")
        adapt.build_python_object(data=self.uniprot_df, early_stopping=early_stopping)
        adapt.write_to_csv_for_admin_import()
