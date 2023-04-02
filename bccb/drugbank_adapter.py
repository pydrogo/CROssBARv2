from __future__ import annotations

from pypath.share import curl, settings, common

from pypath.inputs import drugbank, unichem

from contextlib import ExitStack
from typing import Union
from bioregistry import normalize_curie
from tqdm import tqdm
from time import time
import collections

import pandas as pd
import numpy as np

from biocypher._logger import logger

from enum import Enum, auto

logger.debug(f"Loading module {__name__}.")


class DrugBankNodeField(Enum):
    # primary fields
    DRUGBANK_ID = "drugbank_id"
    NAME = "name"
    CAS_NUMBER = "cas_number"
    GROUPS = "groups"
    GENERAL_REFERENCES = "general_references"
    ATC_CODES = "atc_codes"
    
    # property fields
    INCHI = "InChI"
    INCHIKEY = "InChIKey"
    
    # drugbank external fields
    KEGG_DRUG = "KEGG Drug"
    RXCUI = "RxCUI"
    PHARMGKB = "PharmGKB"
    PDB = "PDB"
    
    # unichem external mappings
    ZINC = "zinc"
    CHEMBL = "chembl"
    BINDINGDB = "bindingdb"
    CLINICALTRIALS = "clinicaltrials"
    CHEBI = "chebi"
    PUBCHEM = "pubchem"
    DRUGCENTRAL = "drugcentral"
    
    @classmethod
    def get_primary_fields(cls):
        return [cls.NAME, cls.CAS_NUMBER, cls.GROUPS, cls.GENERAL_REFERENCES, cls.ATC_CODES, cls.DRUGBANK_ID]
    
    @classmethod
    def get_property_fields(cls):
        return [cls.INCHI, cls.INCHIKEY]
    
    @classmethod
    def get_drugbank_external_fields(cls):
        return [cls.KEGG_DRUG, cls.RXCUI, cls.PHARMGKB, cls.PDB]
    
    @classmethod
    def get_unichem_external_mappings(cls):
        return [cls.ZINC, cls.CHEMBL, cls.BINDINGDB, cls.CLINICALTRIALS, cls.CHEBI, cls.PUBCHEM, cls.DRUGCENTRAL]
    
class PrimaryNodeIdentifier(Enum):
    DRUGBANK = "drugbank"
    KEGG_LIGAND = "kegg_ligand" # kegg ligand id (not drug id)
    DRUGCENTRAL = "drugcentral"
    CHEMBL = "chembl"
    RXCUI = "RxCUI"
    

class DrugBank:
    """
    Class that downloads drug data using pypath and reformats it to be ready
    for import into BioCypher.
    """
    
    def __init__(self, drugbank_user:str, drugbank_passwd:str, node_fields:Union[list[DrugBankNodeField], None] = None,
                primary_node_id:Union[PrimaryNodeIdentifier, None] = None):
        """
        Args
            drugbank_user: drugbank username
            drugbank_passwd: drugbank password
        """
        
        self.user = drugbank_user
        self.passwd = drugbank_passwd
        
        # set node fields
        self.set_node_fields(node_fields=node_fields)
        
        # set primary id of drug nodes
        self.set_primary_id(primary_node_id=primary_node_id)
        
        
    def download_drugbank_data(self, cache=False, debug=False, retries=6,):
        """
        Wrapper function to download drug data from various databases using pypath.

        Args
            cache: if True, it uses the cached version of the data, otherwise
            forces download.
            debug: if True, turns on debug mode in pypath.
            retries: number of retries in case of download error.
        """
        
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())
                
            
            t0 = time()
            
            self.download_drugbank_node_data()
            
            
            t1 = time()
            logger.info(f'All data is downloaded in {round((t1-t0) / 60, 2)} mins'.upper())
            
            
    def download_drugbank_node_data(self) -> None:
        """
        Wrapper function to download DrugBank drug entries, DTI and DDI data using pypath
        Args:
            drugbank_node_fields: node properties that will be included in nodes
        """        
        
        
        logger.debug('Downloading Drugbank drug node data')
        t0 = time()
        
        self.drugbank_data = drugbank.DrugbankFull(user = self.user, passwd = self.passwd)
        
        if self.drugbank_external_node_fields:
            self.drugbank_drugs_external_ids = self.drugbank_data.drugbank_external_ids_full()
        
        if self.property_node_fields:            
            self.drugbank_properties = self.drugbank_data.drugbank_properties_full()
        
        if self.primary_node_fields:
            fields = self.primary_node_fields.copy()
            self.drugbank_drugs_detailed = self.drugbank_data.drugbank_drugs_full(fields = fields)
            
        if self.unichem_mapping_node_fields:
            self.get_unichem_mappings()
        
        t1 = time()
        logger.info(f'Drugbank drug node data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def download_drugbank_dti_data(self):
        
        logger.debug('Downloading Drugbank DTI data')
        t0 = time()
        
        self.drugbank_dti = self.drugbank_data.drugbank_targets_full(fields=['drugbank_id', 'actions', 'references', 'known_action', 'polypeptide',])
        
        t1 = time()
        logger.info(f'Drugbank DTI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def download_drugbank_ddi_data(self):
        logger.debug('Downloading Drugbank DDI data')
        t0 = time()
        
        self.drugbank_ddi = self.drugbank_data.drugbank_drugs_full(fields = "drug_interactions")
        
        t1 = time()
        logger.info(f'Drugbank DDI data is downloaded in {round((t1-t0) / 60, 2)} mins')
        
    def set_node_fields(self, node_fields):
        if node_fields:
            self.all_node_fields = [field.value for field in node_fields]
            self.primary_node_fields = [field.value for field in DrugBankNodeField.get_primary_fields() if field in node_fields]
            self.property_node_fields = [field.value for field in DrugBankNodeField.get_property_fields() if field in node_fields]
            self.drugbank_external_node_fields = [field.value for field in DrugBankNodeField.get_drugbank_external_fields() if field in node_fields]
            self.unichem_mapping_node_fields = [field.value for field in DrugBankNodeField.get_unichem_external_mappings() if field in node_fields]
        else:
            self.all_node_fields = [field.value for field in DrugBankNodeField]
            self.primary_node_fields = [field.value for field in DrugBankNodeField.get_primary_fields()]
            self.property_node_fields = [field.value for field in DrugBankNodeField.get_property_fields()]
            self.drugbank_external_node_fields = [field.value for field in DrugBankNodeField.get_drugbank_external_fields()]
            self.unichem_mapping_node_fields = [field.value for field in DrugBankNodeField.get_unichem_external_mappings()]
        
    def get_unichem_mappings(self) -> None:
        
        self.unichem_mappings = {}
        
        for field in self.unichem_mapping_node_fields: 
            try:
                _dict = unichem.unichem_mapping("drugbank", field)
                
                # remove set data type from value
                _dict = {k:list(v)[0] for k, v in _dict.items()}
                
                # add unichem_mappings
                self.unichem_mappings[field] = _dict
            except TypeError:
                _dict = unichem.unichem_mapping(field, "drugbank")
                
                # switch values and keys
                _dict = {list(v)[0]:k for k, v in _dict.items()}
                
                # add unichem_mappings
                self.unichem_mappings[field] = _dict
                
    def set_primary_id(self, primary_node_id) -> None:
        if primary_node_id:
            self.primary_id = primary_node_id.value
        else:
            self.primary_id = PrimaryNodeIdentifier.DRUGBANK.value
                            
                
    def set_drugbank_to_primary_id_mapping(self) -> None:
        unichem_fields = [field.value for field in DrugBankNodeField.get_unichem_external_mappings()]
        drugbank_external_fields = [field.value for field in DrugBankNodeField.get_drugbank_external_fields()]
        
        if self.primary_id in unichem_fields:
            try:
                self.drugbank_to_primary_id = unichem.unichem_mapping("drugbank", self.primary_id)
                self.drugbank_to_primary_id = {k:list(v)[0] for k, v in drugbank_to_primary_id.items()}
                
            except TypeError:
                self.drugbank_to_primary_id = unichem.unichem_mapping(self.primary_id, "drugbank")
                self.drugbank_to_primary_id = {list(v)[0]:k for k, v in self.drugbank_to_primary_id.items()}
                
        elif self.primary_id in drugbank_external_fields:
            self.drugbank_to_primary_id = {}
            
            for drugbank_id, mapping_dict in self.drugbank_drugs_external_ids.items():
                if mapping_dict.get(self.primary_id, None):
                    self.drugbank_to_primary_id[drugbank_id] = mapping_dict[self.primary_id]
                
                
    def get_drug_nodes(self, label="drug"):
        
        self.node_list = []
                
        if self.primary_id == "drugbank":
            for drug in tqdm(self.drugbank_drugs_detailed):

                node_id = drug.drugbank_id

                props = {}

                if self.primary_node_fields:                
                    props = props | {k:v for k, v in drug._asdict().items() if k in self.primary_node_fields}

                if self.property_node_fields and self.drugbank_properties.get(drug.drugbank_id, None):
                    props = props | {k:v for k, v in self.drugbank_properties[drug.drugbank_id].items() if k in self.property_node_fields}

                if self.drugbank_external_node_fields and self.drugbank_drugs_external_ids.get(drug.drugbank_id, None):
                    props = props | {k:v for k, v in self.drugbank_drugs_external_ids[drug.drugbank_id].items() if k in self.drugbank_external_node_fields}

                if self.unichem_mapping_node_fields:

                    for db, mapping_dict in self.unichem_mappings.items():
                        if mapping_dict.get(drug.drugbank_id, None):
                            props[db] = mapping_dict[drug.drugbank_id]


                self.node_list.append((node_id, label, props))
                
        else:
            # create drugbank to primary id mapping
            self.set_drugbank_to_primary_id_mapping()
            
            for drug in tqdm(self.drugbank_drugs_detailed):
                if self.drugbank_to_primary_id.get(drug.drugbank_id, None):
                    
                    node_id = self.drugbank_to_primary_id[drug.drugbank_id]
                    
                    props = {}
                    
                    if self.primary_node_fields:                
                        props = props | {k:v for k, v in drug._asdict().items() if k in self.primary_node_fields}

                    if self.property_node_fields and self.drugbank_properties.get(drug.drugbank_id, None):
                        props = props | {k:v for k, v in self.drugbank_properties[drug.drugbank_id].items() if k in self.property_node_fields}

                    if self.drugbank_external_node_fields and self.drugbank_drugs_external_ids.get(drug.drugbank_id, None):
                        props = props | {k:v for k, v in self.drugbank_drugs_external_ids[drug.drugbank_id].items() if k in self.drugbank_external_node_fields}

                    if self.unichem_mapping_node_fields:

                        for db, mapping_dict in self.unichem_mappings.items():
                            if mapping_dict.get(drug.drugbank_id, None):
                                props[db] = mapping_dict[drug.drugbank_id]


                    self.node_list.append((node_id, label, props))

            
