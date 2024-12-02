from  ms2pip.single_prediction import SinglePrediction
from pyteomics.parser import parse

def predict_spectrum(peptide : str, charge : int):
    """
    Predicts the spectrum of the peptide string
    """
    peptide = parse(peptide)
    modstring = ""

    for idx, aa in enumerate(peptide):

        match aa:
            case('[carb]C'):
                peptide[idx] = 'C'
                modstring += f'{idx + 1}|Carbamidomethyl|'

            case('[oxid]M'):
                peptide[idx] = 'M'
                modstring += f'{idx + 1}|Oxidation|'

    peptide = ''.join(peptide)

    if modstring == "":
        modstring = "-"
    else:
        modstring = modstring[:-1]

    ms2pip_sp = SinglePrediction(modification_strings=["Carbamidomethyl,57.021464,opt,C", "Oxidation,15.994915,opt,M"])

    mz, intensity, _ = ms2pip_sp.predict(peptide, modstring, charge)

    del ms2pip_sp
    
    return (mz, intensity)
