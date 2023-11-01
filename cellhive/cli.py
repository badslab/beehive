"""Command line interface dispatcher.

"""

import importlib

import click


@click.group()
def cli():
    """Click group."""


@cli.command()
def version():
    """Print version to screen."""
    print(importlib.metadata.version("cellhive"))


def main():
    """Main CLI function."""
    cli()
