# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # Add project root to sys.path

# -- Project information -----------------------------------------------------
project = 'flag'
copyright = '2025, Joe McNease'
author = 'Joe McNease'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_tabs.tabs'
]

templates_path = ['_templates']

# Exclude specific files from scanning entirely
exclude_patterns = [
    "**/co2_tmp.py",
    "**/co2_eos.py",
]

autosummary_generate = True

# -- Autodoc: skip specific modules or members --------------------------------
def skip_members(app, what, name, obj, skip, options):
    """
    Skip autodoc generation for unwanted modules or members.
    """
    # Skip the specific files/modules entirely
    if name in ("flag.co2_tmp", "flag.co2_eos"):
        return True

    # Skip private members (optional)
    if name.startswith("_"):
        return True

    return skip

def setup(app):
    app.connect("autodoc-skip-member", skip_members)

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
