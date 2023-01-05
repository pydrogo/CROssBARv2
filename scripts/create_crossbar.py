"""
CROssBAR generation through BioCypher script
"""

from bccb.uniprot_adapter import (
    Uniprot,
    UniprotNode,
    UniprotNodeField,
    UniprotEdge,
)
from bccb.ppi_adapter import (
    PPI,
    IntactEdgeField,
    BiogridEdgeField,
    StringEdgeField,
)

import biocypher

# Source configuration
uniprot_node_types = [
    UniprotNode.PROTEIN,
    UniprotNode.GENE,
    UniprotNode.ORGANISM,
]

uniprot_node_fields = [
    UniprotNodeField.PROTEIN_SECONDARY_IDS,
    UniprotNodeField.PROTEIN_LENGTH,
    UniprotNodeField.PROTEIN_MASS,
    UniprotNodeField.PROTEIN_ORGANISM,
    UniprotNodeField.PROTEIN_ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_PROTEOME,
    UniprotNodeField.PROTEIN_EC,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS,
    UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS,
    UniprotNodeField.PROTEIN_VIRUS_HOSTS,
    UniprotNodeField.PROTEIN_KEGG_IDS,
]

uniprot_edge_types = [
    UniprotEdge.PROTEIN_TO_ORGANISM,
    UniprotEdge.GENE_TO_PROTEIN,
]

intact_edge_fields = [
    IntactEdgeField.SOURCE,
    IntactEdgeField.PUBMED_IDS,
    IntactEdgeField.INTACT_SCORE,
    IntactEdgeField.METHODS,
    IntactEdgeField.INTERACTION_TYPES,
]

biogrid_edge_fields = [
    BiogridEdgeField.SOURCE,
    BiogridEdgeField.PUBMED_IDS,
    BiogridEdgeField.EXPERIMENTAL_SYSTEM,
]

string_edge_fields = [
    StringEdgeField.SOURCE,
    StringEdgeField.COMBINED_SCORE,
    StringEdgeField.PHYSICAL_COMBINED_SCORE,
]

# Run build
def main():
    """
    Main function. Call individual adapters to download and process data. Build
    via BioCypher from node and edge data.
    """

    driver = biocypher.Driver(
        offline=True,
        db_name="neo4j",
        wipe=True,
        quote_char="'",
        delimiter="\t",
        array_delimiter="|",
        user_schema_config_path="config/schema_config.yaml",
        skip_bad_relationships=True,
        skip_duplicate_nodes=True,
    )

    uniprot_adapter = Uniprot(
        organism="9606",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        test_mode=True,
    )

    uniprot_adapter.download_uniprot_data(
        cache=True,
        retries=5,
    )

    driver.write_nodes(uniprot_adapter.get_nodes())
    driver.write_edges(uniprot_adapter.get_edges())

    # ppi_adapter = PPI(organism=9606,
    #                 intact_fields=intact_edge_fields,
    #                 biogrid_fields=biogrid_edge_fields,
    #                 string_fields=string_edge_fields,
    # )

    # ppi_adapter.download_intact_data()
    # ppi_adapter.intact_process()

    # ppi_adapter.download_biogrid_data()
    # ppi_adapter.biogrid_process()

    # ppi_adapter.download_string_data()
    # ppi_adapter.string_process()

    # ppi_adapter.merge_all()

    # driver.write_edges(ppi_adapter.get_edges())

    driver.write_import_call()
    driver.log_missing_bl_types()
    driver.log_duplicates()
    driver.show_ontology_structure()


if __name__ == "__main__":
    main()
