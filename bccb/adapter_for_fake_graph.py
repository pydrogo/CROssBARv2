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
        db_uri="URİ_HERE",
        db_user="neo4j",
        db_passwd="your_password_here",
        user_schema_config_path="config/schema_config.yaml",
        wipe=False,
        offline=False,
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
        )

        self.nodes = None
        self.edges = None

    def build_python_object(self):
        """
        Loads fake example data.
        """
        nodes = pd.read_csv(r"bccb/data/nodes/nodes.csv")
        edges = pd.read_csv(r"bccb/data/edges/edges.csv")

        self.nodes = nodes
        self.edges = edges

    def translate_python_object_to_neo4j(self, nodes=None, edges=None):
        """
        Loads a pypath network into the biocypher (Neo4j) backend.

        Args:
            - network (pypath.core.network.Network): A network database
              object. If `None`, the value of :py:attr:`network` will be
              used.
        """

        def gen_nodes(nodes):

            for idx, row in tqdm(nodes.iterrows()):
                if isinstance(row["Name"], str):
                    _type = self.process_type(row["Type"])
                    _id = str(row["ID"])
                    _props = {
                        "source_db": str(row["Source Database"]),
                        "name": str(row["Name"]),
                    }
                    yield (_id, _type, _props)

                else:
                    _type = self.process_type(row["Type"])
                    _id = str(row["ID"])
                    _props = {"source_db": str(row["Source Database"])}
                    yield (_id, _type, _props)

        id_type_tuples = gen_nodes(nodes or self.nodes)
        # print(next(id_type_tuples))
        self.bcy.add_nodes(id_type_tuples)

        def gen_edges(edges):

            for idx, row in tqdm(edges.iterrows()):
                _source = str(row["Source ID"])
                _target = str(row["Target ID"])
                _type = str(row["Label"])
                _props = {
                    "source_type": self.process_type(row["Source Type"]),
                    "target_type": self.process_type(row["Target Type"]),
                }
                yield (_source, _target, _type, _props)

        src_tar_type_tuples = gen_edges(edges or self.edges)
        self.bcy.add_edges(src_tar_type_tuples)

    def write_nodes(self, nodes, db_name):
        """
        Writes biocypher nodes to CSV files for admin import.

        Args:
            nodes (list): A list of nodes.
            db_name (str): Name of the database (Neo4j graph) to use.
        """

        # write nodes
        def gen_nodes(nodes):

            for idx, row in tqdm(nodes.iterrows()):
                if isinstance(row["Name"], str):
                    _type = process_type(row["Type"])
                    _id = str(row["ID"])
                    _props = {
                        "source_db": str(row["Source Database"]),
                        "name": str(row["Name"]),
                    }
                    yield (_id, _type, _props)

                else:
                    _type = process_type(row["Type"])
                    _id = str(row["ID"])
                    _props = {"source_db": str(row["Source Database"])}
                    yield (_id, _type, _props)

        id_type_tuples = list(gen_nodes(nodes or self.nodes))

        self.bcy.write_nodes(id_type_tuples, db_name=db_name)

    def write_edges(self, edges, db_name):
        """
        Writes biocypher edges to CSV files for admin import.

        Args:
            edges (list): A list of edges.
            db_name (str): Name of the database (Neo4j graph) to use.
        """

        # write edges
        def gen_edges(edges):

            for idx, row in tqdm(edges.iterrows()):
                _source = str(row["Source ID"])
                _target = str(row["Target ID"])
                _type = str(row["Label"])
                _props = {
                    "source_type": process_type(row["Source Type"]),
                    "target_type": process_type(row["Target Type"]),
                }
                yield (_source, _target, _type, _props)

        src_tar_type_tuples = list(gen_edges(edges or self.edges))

        self.bcy.write_edges(src_tar_type_tuples, db_name=db_name)

    def write_to_csv_for_admin_import(
        self, nodes=None, edges=None, db_name=None
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


def process_type(type_):
    if type_ == "Cell/Tissue":
        type_ = type_.replace("/", "_")
    elif type_ == "Side Effect":
        type_ = type_.replace(" ", "_")

    return str(type_)


def _process_id(identifier):
    """
    Replace critical symbols in ids so that neo4j doesn't throw
    a type error.
    """
    identifier = str(identifier)

    """
    replace_characters = [':', '-']
    
    for character in replace_characters:            
        if character in identifier:
            identifier = identifier.replace(character, "_")
    """

    return identifier
