import typer
from CometSearch.search import main
from typing_extensions import Annotated


app = typer.Typer()

@app.command()
def comet(
    experiment_filename :str,
    protein_database : str, 
    processes : int,
    spectra : int
    ):
    """
    Comet in Python
    """

    main(experiment_filename, protein_database, processes, spectra)


if __name__ == "__main__":
    app()




