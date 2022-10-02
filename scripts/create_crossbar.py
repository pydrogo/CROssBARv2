"""
CROssBAR generation through BioCypher script
"""

# import sys
# sys.path.append('') # fix weird poetry behaviour I don't understand
# may not be necessary depending on your setup

from bccb.protein import Uniprot_data

uniprot_data =  Uniprot_data()
uniprot_data.uniprot_data_download()
uniprot_data.build_dataframe()
uniprot_data.call_biocypher_adapter(early_stopping=5000) # if you want to process whole dataset make early_stopping None
