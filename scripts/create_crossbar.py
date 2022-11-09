"""
CROssBAR generation through BioCypher script
"""

# import sys
# sys.path.append('') # fix weird poetry behaviour I don't understand
# may not be necessary depending on your setup

from bccb.uniprot_adapter import Uniprot

uniprot_data = Uniprot()
uniprot_data.uniprot_data_download(cache=True)
uniprot_data.write_uniprot_nodes_and_edges()

