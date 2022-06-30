"""
CROssBAR generation through BioCypher script
"""

from adapter import BiocypherAdapter

adapt = BiocypherAdapter(wipe=True)
adapt.bcy.current_db

adapt.build_python_object()
adapt.translate_python_object_to_neo4j()
