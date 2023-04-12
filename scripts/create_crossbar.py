"""
CROssBAR generation through BioCypher script
"""

from bccb.uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotEdgeField,
)

from biocypher import BioCypher

# Source configuration
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
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
    UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS,
    UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS,
    UniprotNodeField.PROTEIN_VIRUS_HOSTS,
    UniprotNodeField.PROTEIN_KEGG_IDS,
]

uniprot_edge_types = [
    UniprotEdgeType.PROTEIN_TO_ORGANISM,
    UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_edge_fields = [
    UniprotEdgeField.GENE_ENSEMBL_GENE_ID,
]


# Run build
def main():
    """
    Main function. Call individual adapters to download and process data. Build
    via BioCypher from node and edge data.
    """

    # Start biocypher
    bc = BioCypher()

    # Start uniprot adapter and load data
    uniprot_adapter = Uniprot(
        organism="9606",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        edge_fields=uniprot_edge_fields,
        test_mode=True,
    )

    uniprot_adapter.download_uniprot_data(
        cache=True,
        retries=5,
    )

    # Write uniprot nodes and edges
    bc.write_nodes(uniprot_adapter.get_nodes())

    # Write import call and other post-processing
    bc.write_import_call()
    bc.summary()


if __name__ == "__main__":
    main()
