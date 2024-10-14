import multiprocessing
import time
from pyteomics import mzml, fasta, parser, auxiliary, mass
from CometSearch.utils import fragments, masstocharge_to_dalton, tolerance_bounds
from collections import OrderedDict

pept_index_unsorted = {}
pept_list = []


def add_to_dict(element : str, element_mass : float):
    
    if element_mass in pept_index_unsorted and element not in pept_index_unsorted[element_mass]:
        pept_index_unsorted[element_mass].append(element)
    else:
        pept_index_unsorted[element_mass] = [element]

        

def process_worker(mzml_entry):

    if "MSn spectrum" in mzml_entry:
    #auxiliary.print_tree(mzml_entry)
        for precursor in mzml_entry["precursorList"]["precursor"]:
            if precursor["isolationWindow"]["ms level"] == 1:
                for sel_ion in precursor["selectedIonList"]["selectedIon"]:

                    m_z = sel_ion["selected ion m/z"]
                    charge = sel_ion["charge state"]
                    mass_mzml = masstocharge_to_dalton(m_z, charge)
                    lower, upper = tolerance_bounds(mass_mzml)


def search(filename : str, database : str, processes : int):

    with fasta.read(database) as db:
                    
        for db_entry in db:

            header = db_entry[0]
            sequence = db_entry[1]
            cleaved_sequence = parser.cleave(sequence, "trypsin", 1)

            for element in cleaved_sequence:
                if "X" not in element:
                    element_mass = mass.fast_mass(element)

                    add_to_dict(element, element_mass)

                    if "M" in element:
                        add_to_dict(element, element_mass + 15.994915)

                    if "C" in element:
                        add_to_dict(element, element_mass + 57.021464)

        for peptmass, pepts in pept_index_unsorted.items():
            t = [peptmass]
            for pept in pepts:
                t.append(pept)
            pept_list.append(t)

        #pept_list_sorted = [] = OrderedDict(sorted(pept_index_unsorted.items(), key=lambda x : x[0]))
        pept_list.sort(key=lambda x: x[0])



    with mzml.read(filename) as spectra:
        for protein in spectra:
            process_worker(protein)



""" with mzml.read(filename) as spectra, multiprocessing.Pool(processes) as p:
        #if __name__ == '__main__':
        
            for result in p.imap_unordered(process_worker, spectra):
                print(result) """




                            