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

    def build_python_object(self, data=None):
        """
        Loads uniprot data and preprocess it
        """

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

        id_type_tuples = list(gen_nodes(nodes))

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

        src_tar_type_tuples = list(gen_edges(edges))

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



def gen_nodes(nodes):
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

def gen_edges(edges):
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
