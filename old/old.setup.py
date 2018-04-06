#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup,find_packages
'''

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
'''

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()
#
import pkg_resources

try:
    __version__ = pkg_resources.get_distribution(__name__).version
except:
    __version__ = 'unknown'

setup(
	name = "backmap",
	version = "0.0.1",
	author = "Ranjan Mannige",
	author_email = "ranjanmannige@gmail.com",
	description = ("Module that allows for the visualization of protein backbone structure as graphs"),
	long_description=read('README.md'),
	long_description_content_type='text/markdown',
	license = "MIT",
	keywords = "BackMAP protein backbone dynamics visualization",
	url = "http://www.github.com/ranjanmannige/backmap",
	packages=find_packages(exclude=['biopython','Bio']),
	#packages=['backmap','matplotlib','numpy','copy','string','glob'],
	install_requires=[
	                    'matplotlib',
	                    'numpy'
	                 ], #,'copy','string','glob'
	#tests_requires=['pytest'],
	py_modules=['backmap'],
	test_suite='nose.collector',
    tests_require=['nose'],
	#test_suite="tests",
	classifiers=[
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 3 - Alpha',
		# Indicate who your project is intended for
		'Intended Audience :: Science/Research',
		'Operating System :: OS Independent',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Chemistry',
		'Topic :: Scientific/Engineering :: Visualization'
		# Pick your license as you wish (should match "license" above)
		'License :: OSI Approved :: MIT License',
		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		#'Programming Language :: Python :: 2.7',
		#'Programming Language :: Python :: 3',
		#'Programming Language :: Python :: 3.2',
		#'Programming Language :: Python :: 3.3',
		#'Programming Language :: Python :: 3.4',
		#'Programming Language :: Python :: 3.5',
		#'Programming Language :: Python :: 3.6',
	],
	python_requires='~=2.7', # '~=' is similar to '>=2.7', but without the commitment for the far future (e.g., Python 4)
	#project_urls={
	#'Documentation': 'https://packaging.python.org/tutorials/distributing-packages/',
	#'Funding': 'https://donate.pypi.org',
	#'Say Thanks!': 'http://saythanks.io/to/example',
	#'Source': 'https://github.com/pypa/sampleproject/',
	#'Tracker': 'https://github.com/pypa/sampleproject/issues',
	#},
)
