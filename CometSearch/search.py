import multiprocessing
import time
from pyteomics import mzml, fasta, parser, auxiliary, mass
from CometSearch.utils import fragments, masstocharge_to_dalton, tolerance_bounds

pept_index_unsorted = {}
pept_list = []


def add_to_dict(element : str, element_mass : float):
    
    if element_mass in pept_index_unsorted and element not in pept_index_unsorted[element_mass]:
        pept_index_unsorted[element_mass].append(element)
    else:
        pept_index_unsorted[element_mass] = [element]

def xcorr():
    pass  

def identification(mzml_entry, pep_index):
    bounds = (-1, -1)
    if "MSn spectrum" in mzml_entry:
    #auxiliary.print_tree(mzml_entry)
        for precursor in mzml_entry["precursorList"]["precursor"]:

            if precursor["isolationWindow"]["ms level"] == 1:

                for sel_ion in precursor["selectedIonList"]["selectedIon"]:

                    m_z = sel_ion["selected ion m/z"]
                    charge = sel_ion["charge state"]
                    mass_mzml = masstocharge_to_dalton(m_z, charge)
                    lower, upper = tolerance_bounds(mass_mzml)
                    bounds = (lower, upper)

    return bounds



def main(experiment_filename : str, protein_database : str, processes : int, spectra_amount : int):

    with fasta.read(protein_database) as db:
                    
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

        pept_list.sort(key=lambda x: x[0])


    manager = multiprocessing.Manager()
    pep_index = manager.list(pept_list) # <- peptide index

    total_results = []


    with multiprocessing.Pool(processes) as pool:
        mzml_reader = mzml.read(experiment_filename)

        all_specs_read = False

        while not all_specs_read:
            spec_buffer = []

            for _ in range(spectra_amount): 
                try:
                    spec_buffer.append(next(mzml_reader))
                except StopIteration:
                    all_specs_read = True

            results = [
                pool.apply_async(identification, args=(spec, pep_index))
                for spec in spec_buffer
            ]

            for result in results:
                total_results.append(result.get())
                
    print(total_results)




