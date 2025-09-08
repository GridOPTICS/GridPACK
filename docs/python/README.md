# GridPACK Python Documentation #

## Documentation Generation ##

This is the Sphinx-based documentation for the GridPACK Python
interface.  It is auto-generated, but the `gridpack` module must be
available to Python when it is generated.  It's best to create a
virtual Python environment to provide compatible Sphinx and GridPACK
modules.  At the time of writing, the following Sphinx module are
required:

```
pip install sphinx
pip intall sphinxcontrib-restbuilder
pip install sphinx_rtd_theme
pip install sphinx_subfigure
```

The GridPACK module should also be installed in the environment.  


## Incorporate into Manual Documentation ## 

  * Generate Python interface documentation in RST format:

    ```
    make rst
    ```
  * Copy the generated files in `_build/rst`, except `index.rst`, into
    the user manual documentation tree.  
  * In the user manual document source, modify `python/index.rst` if
    files have been added. Add any new files to the GIT repository.
  * Generate the user manual in HTML format and check.


