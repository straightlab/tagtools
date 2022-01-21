import argparse
import logging

log = logging.getLogger(__name__)

def main(args=None):

    parser = argparse.ArgumentParser(prog='tagtools', description='Command line tool for the creation and quantification of .tag.bam files')
    parser.add_argument(
        '--loglevel', default='info', help='Log level',
        choices=['debug', 'info', 'warning', 'error', 'critical'],
    )
    subparsers = parser.add_subparsers(help='Sub-commands')

    from .tag import add_subcommand_tag
    add_subcommand_tag(subparsers)

    from .resolve import add_subcommand_resolve
    add_subcommand_resolve(subparsers)
    from .quant import add_subcommand_quant
    add_subcommand_quant(subparsers)

    # Parse all command line arguments
    args = parser.parse_args(args)

    # This is not a good way to handle the cases
    # where help should be printed.
    # TODO: there must be a better way?
    if hasattr(args, 'func'):
        # Call the desired subcommand function
        logging.basicConfig(level=args.loglevel.upper())
        args.func(args)
        return 0
    else:
        parser.print_help()
        return 0
