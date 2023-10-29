"""
Behave!
"""

import logging

FORMAT = "%(message)s"
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%X]")

lg = logging.getLogger('termite')

from cellhive.metadata import check, layers, md, obs, obsm
