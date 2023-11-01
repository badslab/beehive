"""
Behave!
"""

import importlib.metadata
import logging

__version__ = importlib.metadata.version("cellhive")


FORMAT = "%(levelname)s %(name)s:%(lineno)s - %(message)s %(asctime)s "
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%X]")

lg = logging.getLogger('termite')
