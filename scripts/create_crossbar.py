"""
CROssBAR generation through BioCypher script
"""

# import sys
# sys.path.append('') # fix weird poetry behaviour I don't understand
# may not be necessary depending on your setup

from bccb.adapter_for_fake_graph import BiocypherAdapter

adapt = BiocypherAdapter(offline=True, db_name="bcv3")

adapt.build_python_object()
adapt.write_to_csv_for_admin_import()
