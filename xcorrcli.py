import typer
from Search.search import main
from typing_extensions import Annotated
from os import cpu_count
from wakepy import keep

app = typer.Typer()

@app.command()
def cli(
    sample_filename : Annotated[str, typer.Argument(help="The name of the sample file")],
    protein_database : Annotated[str, typer.Argument(help="The name of the protein database file")], 
    p : Annotated[int, typer.Option(help="Amount of processes to use in parallel")] = cpu_count() - 2,
    s : Annotated[int, typer.Option(help="Amount of spectra loading at a time")] = 5000,
    ps : Annotated[bool, typer.Option(help="Predict the spectrum")] = False
    ):
    """
    XCorr Reimplementation
    """
    if __name__ == "__main__":
        with keep.presenting():
            main(sample_filename, protein_database, p, s, ps)


if __name__ == "__main__":
    app()




