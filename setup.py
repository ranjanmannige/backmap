#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for project.

    This file was generated with PyScaffold 2.5, a tool that easily
    puts up a scaffold for your new Python project. Learn more under:
    http://pyscaffold.readthedocs.org/

# python 2.7 test
virtualenv bmtest
source bmtest/bin/activate
(to get out: '> deactivate')
python setup.py develop
python setup.py install
python setup.py test
git tag -a v0.1 -m "first"

# python 3.0 test
virtualenv -p python3 bmtest
source bmtest/bin/activate
(to get out: '> deactivate')
python setup.py debug
python setup.py install
python setup.py test
git tag -a v0.1 -m "first"

"""

import sys
from setuptools import setup
import versioneer
#from setuptools_scm import version

# Add here console scripts and other entry points in ini-style format
entry_points = """
[console_scripts]
# script_name = z.module:function
# For example:
# fibonacci = z.skeleton:run
"""

def setup_package():
	needs_sphinx = {'build_sphinx', 'upload_docs'}.intersection(sys.argv)
	
	sphinx = ['sphinx'] if needs_sphinx else []
	setup(setup_requires=[] + sphinx,
                version=versioneer.get_version(),
                cmdclass=versioneer.get_cmdclass(),
			entry_points=entry_points)

if __name__ == "__main__":
    setup_package()
