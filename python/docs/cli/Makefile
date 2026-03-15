# Makefile for Sphinx documentation
#
# Build targets:
#   make html      - Build HTML documentation
#   make latexpdf  - Build PDF documentation via LaTeX
#   make clean     - Remove build artifacts

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help clean html latexpdf

help:
	@echo "Available targets:"
	@echo "  html      - Build HTML documentation"
	@echo "  latexpdf  - Build PDF via LaTeX"
	@echo "  clean     - Remove build artifacts"

clean:
	rm -rf $(BUILDDIR)

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS)
	@echo "HTML docs built in $(BUILDDIR)/html"

latexpdf:
	$(SPHINXBUILD) -b latex $(SOURCEDIR) $(BUILDDIR)/latex $(SPHINXOPTS)
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "PDF built in $(BUILDDIR)/latex/gridpack_cli_manual.pdf"
