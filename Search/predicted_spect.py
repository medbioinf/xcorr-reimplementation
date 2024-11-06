from  ms2pip.single_prediction import SinglePrediction

def predict_spectrum(peptide : str, charge : int):
    ms2pip_sp = SinglePrediction(modification_strings=["Carbamidomethyl,57.021464,opt,C", "Oxidation,15.994915,opt,M"])

    mz, intensity, annotation = ms2pip_sp.predict(peptide, "6|Carbamidomethyl|3|Oxidation", charge)

    return (mz, intensity)


