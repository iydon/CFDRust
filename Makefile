CARGO = cargo
POETRY = poetry

PYTHON = $(POETRY) run python
MATURIN = $(POETRY) run maturin

.PHONY: help init install develop release shell test

help:     ## Print help information
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

init:     ## Initialize the configuration files
	cp config/* .

install:  ## Install Python-side dependencies
	$(POETRY) install

develop:  ## Build the crate and install it as a python module directly in the current virtualenv
	$(MATURIN) develop

release:  ## Build optimized artifacts
	$(MATURIN) develop --release --strip

shell:    ## Spawn a shell within the virtual environment
	$(POETRY) shell

test:     ## Run the tests with cargo
	$(CARGO) test
