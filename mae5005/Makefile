PYTHON = python3

.PHONY: help test

help:   ## Print help information
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

test:   ## Run the unit tests
	$(PYTHON) test.py \
		homework/1/input_and_output.f90 \
		homework/1/bailey_borwein_plouffe_formula.f90

clean:  ## Clear compiled files
	rm *.out
