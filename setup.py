#!/usr/bin/env python

"""
Install superbpp using conda with the command:
    conda install superbpp -c conda-forge -c bioconda

For developers:
    git clone https://github.com/eaton-lab/superbpp
    cd ./superbpp
    conda env create -f environment.yml
    conda activate dev
    pip install -e . --no-deps
"""

import re
from setuptools import setup

# parse version from init.py
with open("velocitree/__init__.py") as init:
    CUR_VERSION = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        init.read(),
        re.M,
    ).group(1)

# setup installation
setup(
    name="velocitree",
    packages=["velocitree"],
    version=CUR_VERSION,
    author="Henry Landis and Deren Eaton",
    author_email="de2356@columbia.edu",
    install_requires=[],
    entry_points={
        'console_scripts': ['velocitree = velocitree.__main__:app']},
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
