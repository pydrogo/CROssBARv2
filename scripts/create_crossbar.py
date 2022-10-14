"""
CROssBAR generation through BioCypher script
"""

# import sys
# sys.path.append('') # fix weird poetry behaviour I don't understand
# may not be necessary depending on your setup

from bccb.protein import Uniprot_data

uniprot_data = Uniprot_data()
uniprot_data.uniprot_data_download(cache=True)
uniprot_data.build_dataframe()
uniprot_data.build_nodes_and_edges(
    early_stopping=50
)  # if you want to process whole dataset make early_stopping None
uniprot_data.call_biocypher_adapter()
