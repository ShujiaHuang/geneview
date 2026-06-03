"""Command-line interface for geneview.

This package provides the CLI layer that allows geneview to be used as a
command-line tool via subcommands (e.g., ``geneview manhattan``, ``geneview qq``).

Each subcommand module exposes two functions:

* ``register(subparsers)`` - add the subcommand's argument parser
* ``run(args)`` - execute the subcommand

Author: Shujia Huang
"""
