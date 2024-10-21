import typer
from Search.search import main
from typing_extensions import Annotated


app = typer.Typer()

@app.command()
def comet(
    sample_filename :str,
    protein_database : str, 
    processes : int,
    spectra : int
    ):
    """
    XCorr Reimplementation
    """
    if __name__ == "__main__":
        main(sample_filename, protein_database, processes, spectra)


if __name__ == "__main__":
    app()




