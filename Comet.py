import typer
from CometSearch.search import search
from typing_extensions import Annotated


app = typer.Typer()

@app.command()
def comet(
    filename :str,
    database : str, #Annotated[str, typer.P(help="to specify a sequence database, overriding entry in parameters file")],
    ):
    """
    Comet in Python
    """

    search(filename, database)


    
if __name__ == "__main__":
    app()




