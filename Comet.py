import typer
from typing_extensions import Annotated

app = typer.Typer()

@app.command()
def comet(
    filename :str,
    p: Annotated[bool, typer.Option(help="to print out a comet.params.new file")] = False,
    q: Annotated[bool, typer.Option(help="to print out a comet.params.new file with more parameter entries")] = False,
    d: Annotated[str, typer.Option(help="to specify a sequence database, overriding entry in parameters file")] = "",
    f: Annotated[int, typer.Option(help="to specify the first/start scan to search, overriding entry in parameters file")] = None,
    l: Annotated[int, typer.Option(help="to specify the last/end scan to search, overriding entry in parameters file")] = None,
    ):
    #print("Usage: 'python Comet.py --help'")
    """
    Comet in Python
    """


    #TODO
    

if __name__ == "__main__":
    app()




