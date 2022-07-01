"""
CROssBAR generation through BioCypher script
"""

import sys
sys.path.append('') # fix weird poetry behaviour I don't understand
# may not be necessary depending on your setup

from bccb.adapter import BiocypherAdapter

adapter = BiocypherAdapter(offline=True)

adapter.build_python_object()
adapter.write_to_csv_for_admin_import()
