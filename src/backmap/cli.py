"""Console script for backmap."""

import typer
import sys
from enum import Enum
from typing_extensions import Annotated
from rich.console import Console

from .backmap import process_PDB, draw_figures

helpme = r"""
============================================
  ____             _    __  __          _____  
 |  _ \           | |  |  \/  |   /\   |  __ \ 
 | |_) | __ _  ___| | _| \  / |  /  \  | |__) |
 |  _ < / _` |/ __| |/ / |\/| | / /\ \ |  ___/ 
 | |_) | (_| | (__|   <| |  | |/ ____ \| |     
 |____/ \__,_|\___|_|\_\_|  |_/_/    \_\_|     
                       (Multi-angle Picture)                                             

This tool provides easily readable "pictures" of protein conformations, 
ensembles, and trajectories saved as either a combined protein databank 
(PDB) structure file, or a directory of such files, and produces graphs.
============================================
MAIN USAGES:
python -m backmap --pdbfn SomeProteinStructure.pdb    [OTHER OPTIONS, IF DESIRED]
python -m backmap --pdbfn /directory/containing/pdbs/ [OTHER OPTIONS, IF DESIRED]
FOR HELP:
python -m backmap --help
============================================"""
if '--help' in sys.argv:
    print(helpme)
app = typer.Typer(add_completion=False, 
                  rich_markup_mode="rich",
                  )
console = Console()

# First option will be the default option
ColortypeOptions = ["Chirality", "SecondaryStructure"]

@app.command(epilog="Made with :heart: by [blue]Ranjan Mannige[/blue]")
def main(pdbfn:Annotated[str,
                         typer.Option("--pdbfn", 
                                     help="""Location of a PDB file or a
                                            directory containing PDB file(s).""")],
         output_dir:Annotated[str,
                              typer.Option("--output-dir",
                              help="""Were the generated figures will be stored.
                              If absent, a report dir will be placed in the
                              PDBs parent directory.""")]='',
         write:Annotated[bool,
                         typer.Option(" /--no-write",
                         help="Don't write figures to the report directory.")] = True,
         show:Annotated[bool, 
                        typer.Option(" /--no-show", 
                        help="Don't show figures while the app is running.")] = True,
         signed: Annotated[bool, 
                           typer.Option("--signed/ ",
                           help="""Use the signed Ramachandran plot.""")] =False,
         colortype:Annotated[str, 
                             typer.Option(
                             case_sensitive=True,
                             help=f"Graph coloring options (options: {[s if ix > 0 else s+' [default]' for ix,s in enumerate(ColortypeOptions)]}).")]
                             =ColortypeOptions[0]
    ):
    #
    if colortype not in ColortypeOptions:
        exit('Must have a valid colortype. The --help lists them.')
    #
    # TODO: identify how to parse multiple files in a single ZIP/TAR archive 
    # (and if it is even useful)
    structure_df = process_PDB(pdbfn=pdbfn,
                               signed=signed)
    #
    draw_figures(structure_df=structure_df, pdbfn=pdbfn, 
                 output_dir=output_dir, write=write, show=show)

if __name__ == "__main__":
    app()
