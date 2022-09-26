
import io
import os
from setuptools import find_packages, setup


VERSION="0.1"


setup(
    name="beehive",
    author="Mark Fiers",
    author_email="mark.fiers@kuleuven.be",
    version=VERSION,
    description="Awesome beehive",
    url="https://github.com/mfiers/beehive/",
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["notebook", "overlay", "secret", "dist",
                                    "tests", ".github"]),
    entry_points={
        "console_scripts": ["beehive = beehive.cli:run"]
    },
)
