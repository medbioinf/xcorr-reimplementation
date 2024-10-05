from pyteomics import mzml, fasta, parser, auxiliary, mass
from CometSearch.utils import fragments, masstocharge_to_dalton

def search(filename : str, database : str):

    with fasta.read(database) as db:
        db_entry = next(db, None)

        while db_entry != None:
            header = db_entry[0]
            sequence = db_entry[1]
            cleaved_sequence = parser.cleave(sequence, "trypsin", 1)

            for element in cleaved_sequence:
                print(element)
                element = element.replace('X', '')
                print(mass.fast_mass(element))

            db_entry = next(db, None)

    with mzml.read(filename) as spectra:
        mzml_entry = next(spectra, None)
        while mzml_entry != None:
            if "MSn spectrum" in mzml_entry:

                #auxiliary.print_tree(mzml_entry)
                m_z = mzml_entry["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
                charge = mzml_entry["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"]

                print(masstocharge_to_dalton(m_z, charge))

            mzml_entry = next(spectra, None)






    