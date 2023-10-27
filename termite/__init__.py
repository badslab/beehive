"""
Behave!
"""

import logging

FORMAT = "%(message)s"
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%X]")

lg = logging.getLogger('termite')

from termite.metadata import layers, md, obsm
