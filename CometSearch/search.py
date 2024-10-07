from pyteomics import mzml, fasta, parser, auxiliary, mass
from CometSearch.utils import fragments, masstocharge_to_dalton, tolerance_check

def search(filename : str, database : str):

    with mzml.read(filename) as spectra:
        
        for mzml_entry in spectra:

            if "MSn spectrum" in mzml_entry:

                #auxiliary.print_tree(mzml_entry)
                m_z = mzml_entry["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
                charge = mzml_entry["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"]
                mass_mzml = masstocharge_to_dalton(m_z, charge)

                print("____________________________")
                print("mzml file entry mass:")
                print(mass_mzml)
                print("Found close masses in fasta:")

                with fasta.read(database) as db:
                    
                    for db_entry in db:

                        header = db_entry[0]
                        sequence = db_entry[1]
                        cleaved_sequence = parser.cleave(sequence, "trypsin", 1)

                        for element in cleaved_sequence:

                    
                            element = element.replace('X', '')
                            db_entry_mass = mass.fast_mass(element)
                    
                            if tolerance_check(mass_mzml, db_entry_mass):

                                print(element)
                                print(db_entry_mass)

                





            



        







    