# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?="-W"
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      ?= _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) --help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

clean:
	rm -rf $(BUILDDIR)/*
	rm -rf auto_examples/

# For local build
local:
	sphinx-build "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	$(SPHINXBUILD) "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

install:
	rm -rf _build/doctrees _build/nimgen.github.io
	# first clone the gh-pages branch because it may ask
	# for password and we don't want to delay this long build in
	# the middle of it
	# --no-checkout just fetches the root folder without content
	# --depth 1 is a speed optimization since we don't need the
	# history prior to the last commit
	git clone -b gh-pages --no-checkout --depth 1 git@github.com:juaml/nimgen.git _build/nimgen.github.io
	touch _build/nimgen.github.io/.nojekyll
	make local
	cd _build/ && \
	cp -r html/* nimgen.github.io && \
	cd nimgen.github.io && \
	git add * && \
	git add .nojekyll && \
	git commit -a -m 'Make install' && \
	git push