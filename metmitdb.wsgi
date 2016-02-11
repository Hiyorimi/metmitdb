activate_this = '/home/metmitdb/metmitdb_venv/bin/activate_this.py'
execfile(activate_this, dict(__file__=activate_this))
import sys
import logging
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, '/home/metmitdb/metmitdb_venv')
sys.path.insert(0, '/home/metmitdb/metmitdb/')
from metmitdb import app as application
