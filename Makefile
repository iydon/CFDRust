CARGO = cargo
POETRY = poetry

PYTHON = $(POETRY) run python
MATURIN = $(POETRY) run maturin

.PHONY: help install develop shell

help:     ## Print help information
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

install:  ## Install Python-side dependencies
	$(POETRY) install

develop:  ## Build the crate and install it as a python module directly in the current virtualenv
	$(MATURIN) develop

shell:    ## Spawn a shell within the virtual environment
	$(POETRY) shell
