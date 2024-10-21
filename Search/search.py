import multiprocessing
import time
from pyteomics import mzml, fasta, parser, auxiliary, mass
from Search.utils import fragments, masstocharge_to_dalton, tolerance_bounds, binary_search


pept_index_unsorted = {}
pept_list = []


def add_to_dict(element : str, element_mass : float):
    
    if element_mass in pept_index_unsorted and element not in pept_index_unsorted[element_mass]:
        pept_index_unsorted[element_mass].append(element)
    else:
        pept_index_unsorted[element_mass] = [element]



def identification(mzml_entry, pep_index, list_length):

    if "MSn spectrum" in mzml_entry:
    #auxiliary.print_tree(mzml_entry)
        for precursor in mzml_entry["precursorList"]["precursor"]:

            if precursor["isolationWindow"]["ms level"] == 1:

                for sel_ion in precursor["selectedIonList"]["selectedIon"]:

                    m_z = sel_ion["selected ion m/z"]
                    charge = sel_ion["charge state"]
                    mass_mzml = masstocharge_to_dalton(m_z, charge)
                    lower, upper = tolerance_bounds(mass_mzml)

                    lower_index = binary_search(pep_index, lower, list_length)

                    if lower_index == -1:
                        return None

                    upper_index = binary_search(pep_index, upper, list_length)

                    if upper_index == -1:
                        upper_index = list_length
                    else:
                        upper_index = upper_index - 1

                    pep_index_slice = pep_index[lower_index:upper_index]

                    return pep_index_slice
    return None



def main(sample_filename : str, protein_database : str, processes : int, spectra_amount : int):
    start = time.time()
    with fasta.read(protein_database) as db:
                    
        for db_entry in db:

            header = db_entry[0]
            sequence = db_entry[1]
            cleaved_sequence = parser.cleave(sequence, "trypsin", 1)

            for element in cleaved_sequence:

                if "X" not in element:
                    element_mass = mass.fast_mass(element)

                    add_to_dict(element, element_mass)
                    m_count = element.count("M")
                    
                    if m_count != 0:

                        for cnt in range(min(m_count, 3)):
                            add_to_dict(element, element_mass + (cnt * 15.994915)) #max 3 modified M

                    c_count = element.count("C")

                    if c_count != 0:
            
                        add_to_dict(element, element_mass + (c_count * 57.021464)) #all C modified

        for peptmass, pepts in pept_index_unsorted.items():
            t = [peptmass]
            for pept in pepts:
                t.append(pept)
            pept_list.append(t)

        pept_list.sort(key=lambda x: x[0])
        print("Peptide Database List created!")
        
    list_length = len(pept_list)

    manager = multiprocessing.Manager()
    pep_index = manager.list(pept_list) # <- peptide index

    total_results = []


    with multiprocessing.Pool(processes) as pool, open("output.txt", "w") as outfile:
        mzml_reader = mzml.read(sample_filename)

        all_specs_read = False

        while not all_specs_read:
            spec_buffer = []

            for _ in range(spectra_amount): 
                try:
                    spec_buffer.append(next(mzml_reader))
                except StopIteration:
                    all_specs_read = True

            results = [
                pool.apply_async(identification, args=(spec, pep_index, list_length))
                for spec in spec_buffer
            ]

            for result in results:
                res = result.get()
                if res is not None:
                    print(res, file=outfile)
                    #total_results.append(res)
    end = time.time()
    print(end - start)




