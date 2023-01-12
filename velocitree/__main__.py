#!/usr/bin/env python

"""
Command line interface for velocitree

velocitree generate --tree ... --model ... ...

velocitree fit --tree ... --model 'all'

"""

from pathlib import Path
import typer
from velocitree import __version__
from velocitree.speciation.regression import Models
from velocitree.speciation.generative import RandomTree, UserTree
from velocitree.speciation.viz import heatmap_tree_plot


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
app = typer.Typer(add_completion=False, context_settings=CONTEXT_SETTINGS)


def version_callback(value: bool):
    "Adding a --version option to the CLI"
    if value:
        typer.echo(f"velocitree {__version__}")
        raise typer.Exit()

@app.callback()
def main(
    version: bool = typer.Option(
        None, 
        "--version", "-v", 
        callback=version_callback, 
        is_eager=True, 
        help="print version and exit."),
    ):
    """
    Call kmerkit commands to access tools in the kmerkit toolkit,
    and kmerkit COMMAND -h to see help options for each tool
    (e.g., kmerkit count -h)
    """
    typer.secho(
        f"velocitree (v.{__version__}): quantify RI on trees",
        fg=typer.colors.MAGENTA, bold=version,
    )

@app.command()
def draw(
    tree: Path = typer.Option(..., help="Path to a newick tree file"),
    ridata: Path = typer.Option(100),
    spdata: Path = typer.Option(None),
    format: str = typer.Option(None),
    ):
    """Generates plots of observed hybrids/RI given a phylogeny.
    """
    typer.secho("Generating plots", fg=typer.colors.MAGENTA)
    dist_RI_scatterplot(subsample)
    heatmap_tree_plot(data.tree, data.clades, subsample)


@app.command()
def generate(
    tree: str = typer.Option(..., help="Path to a newick tree file"),
    ncrosses: int = typer.Option(100),
    ntips: int = typer.Option(100, ),
    nclades: int = typer.Option(1, help="N hierarchical groups in tree"),
    model: Models = typer.Option(Models.linear, help="regression model"),
    name_prefix: str = typer.Option('./test', help="Path prefix for output files"),
    seed: int = typer.Option(12345, help="random seed")
    ):
    """
    Generate a test crossing dataset on an input tree or on a random 
    birth-death tree for performing power analyses.
    """
    if tree:
        typer.secho(
            "Generating a test data set on user input tree", 
            fg=typer.colors.MAGENTA,
        )

    else:
        typer.secho(
            "Generating {ncrosses} observations on a {ntips} tip b-d tree", 
            fg=typer.colors.MAGENTA,
        )
        tool = RandomTree(ntips, nclades, model, seed)
        for param in tool.params:
            print("{}: {}".format(param, tool.params[param]))
        print(tool.spdata.head())
        print(tool.sample_observations(ncrosses).head())
