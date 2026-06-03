"""Main entry point for the geneview command-line interface.

Usage::

    $ geneview <subcommand> [options]

Run ``geneview --help`` for a list of available subcommands, or
``geneview <subcommand> --help`` for subcommand-specific options.

Author: Shujia Huang
"""
import sys
import argparse


__all__ = ["main"]


def _get_version():
    """Retrieve the geneview package version string."""
    try:
        from geneview import __version__
        return __version__
    except (ImportError, AttributeError):
        return "unknown"


def main(argv=None):
    """Main CLI entry point.

    Parameters
    ----------
    argv : list of str or None, optional
        Command-line arguments. If None (default), ``sys.argv[1:]`` is used.
        This parameter is mainly useful for programmatic invocation and testing.

    Returns
    -------
    int
        Exit code: 0 on success, non-zero on failure.
    """
    parser = argparse.ArgumentParser(
        prog="geneview",
        description="Geneview: A genomics data visualization toolkit. "
                    "Provides convenient commands for creating publication-quality "
                    "genomics figures from the command line.",
        epilog="Use 'geneview <command> -h' for more information about a specific command.",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s " + _get_version(),
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="Available visualization commands",
        help="Run 'geneview <command> -h' for command-specific help.",
    )

    # Lazy-import subcommand modules and register their parsers.
    from . import manhattan
    from . import qq
    from . import venn
    from . import admixture
    from . import tracks

    manhattan.register(subparsers)
    qq.register(subparsers)
    venn.register(subparsers)
    admixture.register(subparsers)
    tracks.register(subparsers)

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    try:
        args.func(args)
    except Exception as e:
        sys.stderr.write("[ERROR] %s\n" % str(e))
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
