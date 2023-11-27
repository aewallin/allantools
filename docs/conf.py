# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AllanTools'
copyright = '2023, Anders E. Wallin'
author = 'Anders E. Wallin'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# 'sphinx.ext.napoleon',
extensions = ['sphinx.ext.autodoc',
                'sphinx_rtd_theme']

templates_path = ['_templates']
exclude_patterns = []

root_doc = 'index'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'classic'
# https://sphinx-rtd-theme.readthedocs.io/en/stable/installing.html
html_theme = "sphinx_rtd_theme"
# html_static_path = ['_static']
