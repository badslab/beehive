
import os

import yaml
from pathlib import Path


BASEDIR = Path(os.environ['BEEHIVE_BASEDIR'])
# BASEDIR = Path(__file__).parent.parent.parent
BOKEHDIR = BASEDIR / 'bokeh'
DATADIR = BASEDIR / 'data'
