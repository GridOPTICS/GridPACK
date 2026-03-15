# Configuration file for the Sphinx documentation builder.
#
# GridPACK Command-Line Interface — User Manual

project = 'GridPACK Command-Line Interface'
copyright = '2026, Battelle Memorial Institute'
author = 'GridPACK Dev Team'
release = '3.5'
version = '3.5'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = ['_build']

# HTML output
html_theme = 'alabaster'
html_static_path = ['_static']
html_title = 'GridPACK CLI %s' % release

# LaTeX/PDF output
latex_engine = 'pdflatex'
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '11pt',
    'preamble': r'''
\usepackage{enumitem}
\setlistdepth{99}
''',
}

latex_documents = [
    ('index', 'gridpack_cli_manual.tex',
     'GridPACK Command-Line Interface',
     'GridPACK Dev Team', 'manual'),
]
