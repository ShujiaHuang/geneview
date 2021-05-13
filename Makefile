export SHELL := /bin/bash

test:
	pytest --doctest-modules --cov=geneview --cov-config=.coveragerc geneview

unittests:
	pytest --cov=geneview --cov-config=.coveragerc geneview

lint:
	flake8 geneview
