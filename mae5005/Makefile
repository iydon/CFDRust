PYTHON = python3

.PHONY: help test

help:   ## Print help information
	@fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

test:   ## Run the unit tests
	$(PYTHON) test.py \
		code/1/input_and_output.f90 \
		code/1/bailey_borwein_plouffe_formula.f90 \
		code/2/unidirectional_flow.f90 \
		code/3/steady_state_heat_conduction.f90 \
		code/4/1d_transient_laminar_flow.f90

clean:  ## Clear compiled files
	rm *.dat *.mod *.out
