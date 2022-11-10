"""
CROssBAR generation through BioCypher script
"""

from bccb.uniprot_adapter import Uniprot, Fields
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


fields = [
    Fields.PROTEIN_LENGTH,
    Fields.PROTEIN_MASS,
    Fields.PROTEIN_ORGANISM,
]
uniprot_adapter = Uniprot(
    organism="9606",
    fields=fields,
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
