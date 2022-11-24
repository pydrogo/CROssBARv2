"""
CROssBAR generation through BioCypher script
"""

from bccb.uniprot_adapter import Uniprot, UniprotNodeFields, UniprotEdgeFields
import biocypher

driver = biocypher.Driver(
            offline=True,
            db_name="neo4j",
            wipe=True,
            quote_char="'",
            delimiter="\t",
            array_delimiter="|",
            user_schema_config_path="config/schema_config.yaml",
            skip_bad_relationships=True,
            skip_duplicate_nodes=True
        )


node_fields = [
    UniprotNodeFields.PROTEIN_SECONDARY_IDS,
    UniprotNodeFields.PROTEIN_LENGTH,
    UniprotNodeFields.PROTEIN_MASS,
    UniprotNodeFields.PROTEIN_ORGANISM,
    UniprotNodeFields.PROTEIN_ORGANISM_ID,
    UniprotNodeFields.PROTEIN_NAMES,
    UniprotNodeFields.PROTEOME,
    UniprotNodeFields.EC,
    UniprotNodeFields.GENE_NAMES,
    UniprotNodeFields.ENSEMBL_GENE_IDS,
    UniprotNodeFields.ENTREZ_GENE_IDS,
    UniprotNodeFields.VIRUS_HOSTS,
    UniprotNodeFields.KEGG_IDS,
]

edge_fields = [
    UniprotEdgeFields.PROTEIN_TO_ORGANISM,
    UniprotEdgeFields.GENE_TO_PROTEIN,
]

uniprot_adapter = Uniprot(
    organism="9606",
    node_fields=node_fields,
    edge_fields=edge_fields,
)

uniprot_adapter.download_uniprot_data(
    cache=True, 
    retries=5,
)

driver.write_nodes(uniprot_adapter.get_uniprot_nodes())
driver.write_edges(uniprot_adapter.get_uniprot_edges())

driver.write_import_call()
driver.log_missing_bl_types()
driver.log_duplicates()
driver.show_ontology_structure()
