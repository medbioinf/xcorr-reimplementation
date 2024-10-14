import typer
from CometSearch.search import search
from typing_extensions import Annotated


app = typer.Typer()

@app.command()
def comet(
    filename :str,
    database : str, 
    processes : int
    ):
    """
    Comet in Python
    """

    search(filename, database, processes)


if __name__ == "__main__":
    app()




