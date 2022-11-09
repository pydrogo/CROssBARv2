"""
CROssBAR generation through BioCypher script
"""

from bccb.uniprot_adapter import Uniprot

uniprot_adapter = Uniprot()
uniprot_adapter.uniprot_data_download(cache=True)
uniprot_adapter.write_uniprot_nodes_and_edges()

