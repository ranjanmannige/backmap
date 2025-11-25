# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BackMAP'
copyright = '2025, Ranjan Mannige'
author = 'Ranjan Mannige'
release = '1.0.0'

import os
import sys
#sys.path.insert(0, os.path.abspath('../'))
#sys.path.insert(0, os.path.abspath('../')) 
sys.path.insert(0, os.path.abspath('../src/'))
#sys.path.insert(0, os.path.abspath('../src/backmap/')) 

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx_markdown_builder',
    'sphinx_rtd_theme',
    'myst_parser',
    'nbsphinx',
    'nbsphinx_link'
]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_member_order = 'bysource'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
#html_theme = 'alabaster'

html_static_path = ['_static']
