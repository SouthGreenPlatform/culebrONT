#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import culebrONT
from culebrONT.snakeWrapper.global_variable import *

# The short X.Y version.
version = culebrONT.__version__
# The full version, including alpha/beta/rc tags
release = culebrONT.__version__

rst_prolog = f"""
.. |tools_path| replace:: {GIT_TOOLS_PATH}
"""

# -- Project information -----------------------------------------------------
# General information about the project.
project = 'CulebrONT'
copyright = '2019-2022, J Orjuela (IRD), A Comte (IRD), S Ravel (CIRAD), F Charriat (INRAE), T Vi (IRD, AGI), F Sabot (IRD) and S Cunnac (IRD)'
github_doc_root = 'https://github.com/SouthGreenPlatform/CulebrONT_pipeline/tree/master/docs/'
issues_github_path = 'https://github.com/SouthGreenPlatform/CulebrONT_pipeline/issues'

latex_authors = '''
Julie Orjuela (IRD),\\\\
Aurore Comte (IRD),\\\\
Sebastien Ravel (CIRAD),\\\\
Florian Charriat (INRAE)\\\\
Tram Vi(IRD, AGI),\\\\
Francois Sabot(IRD),\\\\
Sébastien Cunnac(IRD)
'''

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
              'sphinx.ext.viewcode',
              'sphinx.ext.autosectionlabel',
              'sphinx_copybutton',
              'sphinx_rtd_theme',
              'sphinx_click'
              ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = ['.rst', "md"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

master_doc = 'index'

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'analytics_id': 'UA-172723859-1',  #  Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    # 'style_nav_header_background': 'cyan',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': False,
    'titles_only': False
}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "CulebrONT"

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "CulebrONT"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_images/culebront_logo.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_images/culebront_logo2.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If false, no index is generated.
html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True


# -- Options for LaTeX output ---------------------------------------------
# latex_engine = 'pdflatex'

# latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'a4paper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '12pt',
    # Latex figure (float) alignment
    # 'figure_align':'htbp',
    # 'extraclassoptions': 'openany',
    # 'preamble': r'''
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # %%%add number to subsubsection 2=subsection, 3=subsubsection
        # \setcounter{secnumdepth}{0}
        # %%%% Table of content upto 2=subsection, 3=subsubsection
        # \setcounter{tocdepth}{2}
    # ''',

    # 'sphinxsetup': \
        # 'hmargin={0.7in,0.7in}, vmargin={0.7in,0.7in}, \
        # marginpar=1in, \
        # verbatimwithframe=False, \
        # TitleColor={rgb}{0,0,0}, \
        # HeaderFamily=\\rmfamily\\bfseries, \
        # InnerLinkColor={rgb}{0,0,1}, \
        # OuterLinkColor={rgb}{0,0,1}',
# }

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).

# latex_documents = [
  # ('index', 'CulebrONT.tex', 'Documentation',
   # latex_authors, 'manual', True),
# ]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = '_images/culebront_logo.png'

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True

# latex_toplevel_sectioning = 'section'
