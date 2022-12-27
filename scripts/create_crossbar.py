"""
CROssBAR generation through BioCypher script
"""

from bccb.uniprot_adapter import Uniprot, UniprotNodeFields, UniprotEdgeFields
from bccb.ppi_adapter import PPI, IntactEdgeFields, BiogridEdgeFields, StringEdgeFields

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


uniprot_node_fields = [
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

uniprot_edge_fields = [
    UniprotEdgeFields.PROTEIN_TO_ORGANISM,
    UniprotEdgeFields.GENE_TO_PROTEIN,
]

intact_edge_fields = [IntactEdgeFields.SOURCE,
                      IntactEdgeFields.UNIPROT_A,
                      IntactEdgeFields.UNIPROT_B,
                      IntactEdgeFields.PUBMED_IDS,
                      IntactEdgeFields.INTACT_SCORE,
                      IntactEdgeFields.METHODS,
                      IntactEdgeFields.INTERACTION_TYPES,
                     ]
biogrid_edge_fields = [BiogridEdgeFields.SOURCE,
                      BiogridEdgeFields.UNIPROT_A,
                      BiogridEdgeFields.UNIPROT_B,
                      BiogridEdgeFields.PUBMED_IDS,
                      BiogridEdgeFields.EXPERIMENTAL_SYSTEM,
                      ]

string_edge_fields = [StringEdgeFields.SOURCE,
                      StringEdgeFields.UNIPROT_A,
                      StringEdgeFields.UNIPROT_B,
                      StringEdgeFields.COMBINED_SCORE,
                      StringEdgeFields.PHYSICAL_COMBINED_SCORE,
                     ]


uniprot_adapter = Uniprot(
    organism="9606",
    node_fields=uniprot_node_fields,
    edge_fields=uniprot_edge_fields,
)

uniprot_adapter.download_uniprot_data(
    cache=True, 
    retries=5,
)

ppi_adapter = PPI(organism=9606, 
                  intact_fields=intact_edge_fields,
                  biogrid_fields=biogrid_edge_fields, 
                  string_fields=string_edge_fields,
)

ppi_adapter.download_intact_data()
ppi_adapter.intact_process()

ppi_adapter.download_biogrid_data()
ppi_adapter.biogrid_process()

ppi_adapter.download_string_data()
ppi_adapter.string_process()

ppi_adapter.merge_all()

driver.write_nodes(uniprot_adapter.get_uniprot_nodes())
driver.write_edges(uniprot_adapter.get_uniprot_edges())
driver.write_edges(ppi_adapter.get_ppi_edges())

driver.write_import_call()
driver.log_missing_bl_types()
driver.log_duplicates()
driver.show_ontology_structure()
