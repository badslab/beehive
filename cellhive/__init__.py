"""
Behave!
"""

import logging

FORMAT = "%(levelname)s %(name)s:%(lineno)s - %(message)s %(asctime)s "
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%X]")

lg = logging.getLogger('termite')

from cellhive.metadata import check, layers, md, obs, obsm
