from  ms2pip.single_prediction import SinglePrediction

def predict_spectrum(peptide : str, charge : int):
    """
    Predicts the spectrum of the peptide string
    """

    modstring = ""

    for idx, aa in enumerate(peptide):

        match aa:
            case('C'):
                modstring += f'{idx + 1}|Carbamidomethyl|'

            case('M'):
                modstring += f'{idx + 1}|Oxidation|'

    if modstring == "":
        modstring = "-"
    else:
        modstring = modstring[:-1]

    ms2pip_sp = SinglePrediction(modification_strings=["Carbamidomethyl,57.021464,opt,C", "Oxidation,15.994915,opt,M"])

    mz, intensity, annotation = ms2pip_sp.predict(peptide, modstring, charge)

    return (mz, intensity)
