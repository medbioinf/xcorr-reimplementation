from pyteomics import mzml, fasta, parser, auxiliary, mass
from CometSearch.utils import fragments, masstocharge_to_dalton, tolerance_check

def search(filename : str, database : str):

    with mzml.read(filename) as spectra:
        
        for mzml_entry in spectra:

            if "MSn spectrum" in mzml_entry:
                #auxiliary.print_tree(mzml_entry)
                for precursor in mzml_entry["precursorList"]["precursor"]:
                    if precursor["isolationWindow"]["ms level"] == 1:
                        for sel_ion in precursor["selectedIonList"]["selectedIon"]:

                            m_z = sel_ion["selected ion m/z"]
                            charge = sel_ion["charge state"]
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
                                        if "X" not in element:

                                            temp = []
                                            temp.append(mass.fast_mass(element))

                                            if "M" in element:
                                                temp.append(mass.fast_mass(element) + 15.994915)

                                            if "C" in element:
                                                temp.append(mass.fast_mass(element) + 57.021464)

                                            for db_entry_mass in temp:
                    
                                                if tolerance_check(mass_mzml, db_entry_mass):

                                                    print(element)
                                                    print(db_entry_mass) 