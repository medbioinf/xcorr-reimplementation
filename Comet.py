import typer
from CometSearch.search import search
from typing_extensions import Annotated


app = typer.Typer()

@app.command()
def comet(
    filename :str,
    database : str, #Annotated[str, typer.P(help="to specify a sequence database, overriding entry in parameters file")],
    p: Annotated[bool, typer.Option(help="to print out a comet.params.new file")] = False,
    q: Annotated[bool, typer.Option(help="to print out a comet.params.new file with more parameter entries")] = False,
    ap: Annotated[str, typer.Option(help="to specify an alternate parameters file (default comet.params)")] = "",
    n: Annotated[str, typer.Option(help="to specify an alternate output base name; valid only with one input file")] = "",
    f: Annotated[int, typer.Option(help="to specify the first/start scan to search, overriding entry in parameters file")] = None,
    l: Annotated[int, typer.Option(help="to specify the last/end scan to search, overriding entry in parameters file")] = None,
    i: Annotated[str, typer.Option(help="create peptide index file only (specify .idx file as database for index search)")] = ""
    ):
    """
    Comet in Python
    """

    search(filename, database)

    if p:
        pass

    if q:
        pass

    if ap != "":
        pass

    if n != "":
        pass

    if f != None:
        pass

    if l != None:
        pass

    if i != "":
        pass

    
if __name__ == "__main__":
    app()




