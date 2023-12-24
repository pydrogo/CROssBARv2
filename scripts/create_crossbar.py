from bccb.uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotEdgeField,
)

from bccb.ppi_adapter import (
    PPI
)

from bccb.interpro_adapter import (
    InterPro
)

from bccb.go_adapter import (
    GO
)

from bccb.drug_adapter import (
    Drug
)

from bccb.compound_adapter import (
    Compound
)

from bccb.orthology_adapter import (
    Orthology
)

from bccb.disease_adapter import (
    Disease
)

from bccb.phenotype_adapter import (
    HPO
)

from bccb.pathway_adapter import (
    Pathway
)

from biocypher import BioCypher

bc = BioCypher(biocypher_config_path= r"config/biocypher_config.yaml",
               schema_config_path= r"config/schema_config.yaml",
)

# Whether to cache data by pypath for future usage
CACHE = True

# Flag for exporting node and edge files as csv format
export_as_csv = True

# dirs
output_dir_path = "YOUR_PATH"

# user and passwd
drugbank_user = "YOUR_DRUGBANK_USER"
drugbank_passwd = "YOUR_DRUGBANK_PASSWD"

# uniprot configuration
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
]

uniprot_node_fields = [
    UniprotNodeField.PRIMARY_GENE_NAME,
    UniprotNodeField.PROTEIN_LENGTH,
    UniprotNodeField.PROTEIN_MASS,
    UniprotNodeField.PROTEIN_ORGANISM,
    UniprotNodeField.PROTEIN_ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_PROTEOME,
    UniprotNodeField.PROTEIN_SUBCELLULAR_LOCATION,
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
    UniprotEdgeField.GENE_ENTREZ_ID,
]

uniprot_adapter = Uniprot(
        organism="*",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        edge_fields=uniprot_edge_fields,
        test_mode=False,
    )

uniprot_adapter.download_uniprot_data(cache=CACHE, retries=6)

uniprot_nodes = uniprot_adapter.get_nodes()
uniprot_edges = uniprot_adapter.get_edges()

bc.write_nodes(uniprot_nodes)
bc.write_edges(uniprot_edges)
bc.write_import_call()

if export_as_csv:
    uniprot_adapter.export_data_to_csv(path=output_dir_path,
                                    node_data=uniprot_nodes,
                                    edge_data=uniprot_edges)

# PPI
ppi_adapter = PPI(organism=None, 
                  cache=CACHE, 
                  output_dir=output_dir_path,
                  export_csv=export_as_csv)

ppi_adapter.download_intact_data()
ppi_adapter.download_biogrid_data()
ppi_adapter.download_string_data()

ppi_adapter.intact_process()
ppi_adapter.biogrid_process()
ppi_adapter.string_process()

bc.write_edges(ppi_adapter.get_ppi_edges())

# protein domain
interpro_adapter = InterPro(cache=CACHE)

interpro_adapter.download_domain_node_data()
interpro_adapter.download_domain_edge_data()

if export_as_csv:
    interpro_adapter.export_as_csv(path=output_dir_path, 
                                node_csv_name="domain",
                                edge_csv_name="protein_has_domain")

bc.write_nodes(interpro_adapter.get_interpro_nodes())
bc.write_edges(interpro_adapter.get_interpro_edges())

# gene ontology
go_adapter = GO(organism="*")
go_adapter.download_go_data(cache=CACHE)
bc.write_nodes(go_adapter.get_go_nodes())
bc.write_edges(go_adapter.get_go_edges())
if export_as_csv:
    go_adapter.export_as_csv(path=output_dir_path)

# drug
drug_adapter = Drug(drugbank_user=drugbank_user, drugbank_passwd=drugbank_passwd,
                    export_csv=export_as_csv, output_dir=output_dir_path)
drug_adapter.download_drug_data(cache=CACHE)
drug_adapter.process_drug_data()
bc.write_nodes(drug_adapter.get_drug_nodes())
bc.write_edges(drug_adapter.get_dti_edges())
bc.write_edges(drug_adapter.get_ddi_edges())
bc.write_edges(drug_adapter.get_dgi_edges())

# compound
compound_adapter = Compound(export_csv=export_as_csv, output_dir=output_dir_path)
compound_adapter.download_compound_data(cache=CACHE)
compound_adapter.process_compound_data()
bc.write_nodes(compound_adapter.get_compound_nodes())
bc.write_edges(compound_adapter.get_cti_edges())

# orthology
orthology_adapter = Orthology(export_csv=export_as_csv, output_dir=output_dir_path)
orthology_adapter.download_orthology_data(cache=CACHE)
bc.write_edges(orthology_adapter.get_orthology_edges())

# disease
disease_adapter = Disease(drugbank_user=drugbank_user, drugbank_passwd=drugbank_passwd,
                          export_csv=export_as_csv, output_dir=output_dir_path)
disease_adapter.download_disease_data(cache=CACHE)
bc.write_nodes(disease_adapter.get_nodes())
bc.write_edges(disease_adapter.get_edges())

# phenotype
phenotype_adapter = HPO(export_csv=export_as_csv, output_dir=output_dir_path)
phenotype_adapter.download_hpo_data(cache=CACHE)
bc.write_nodes(phenotype_adapter.get_nodes())
bc.write_edges(phenotype_adapter.get_edges())

# pathway
pathway_adapter = Pathway(drugbank_user=drugbank_user, drugbank_passwd=drugbank_passwd)
pathway_adapter.download_pathway_data(cache=CACHE)
pathway_adapter.get_nodes()
pathway_adapter.get_edges()

# Write import call and other post-processing
bc.write_import_call()
bc.summary()
