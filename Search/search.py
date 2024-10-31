# std import
from collections import defaultdict
import json
import multiprocessing
from pathlib import Path
import time
from typing import DefaultDict, Set, List, TextIO, Tuple

# 3rd party import
from pyteomics import mzml, parser, auxiliary, mass
from pyteomics.fasta import read as read_fasta
import numpy as np

# internal imports
from Search.utils import fragments, masstocharge_to_dalton, tolerance_bounds, binary_search
from Search.binning import binning

PEPTIDE_MIN_LENGTH = 6
PEPTIDE_MAX_LENGTH = 50
MAX_MISSED_CLEAVAGES = 2

def create_pept_index(fasta_content: TextIO) -> List[Tuple[float, Tuple[str]]]:
    pept_index: DefaultDict[float, Set[str]] = defaultdict(set)

    with read_fasta(fasta_content) as fasta:
        for (_, sequence) in fasta:
            peptides = parser.cleave(sequence, "trypsin", missed_cleavages=MAX_MISSED_CLEAVAGES, min_length=PEPTIDE_MIN_LENGTH, max_length=PEPTIDE_MAX_LENGTH)
            for peptide in peptides:
                if "X" not in peptide:
                    pep_mass = mass.fast_mass(peptide)

                    c_count = peptide.count("C")
                    static_modified_pep_mass = pep_mass + (c_count * 57.021464)

                    pept_index[static_modified_pep_mass].add(peptide) # add static modified sequence
                    
                    m_count = peptide.count("M")
                    
                    if m_count != 0:
                        for cnt in range(min(m_count, 3)):
                            pept_index[static_modified_pep_mass + (cnt * 15.994915)].add(peptide) # add modified M

    print("Read done!")
    pept_index: List[Tuple[float, Tuple[str]]] = [
        (mass, tuple(peptides))
        for mass, peptides in pept_index.items()
    ]

    print("Conversion done!")

    pept_index.sort(key=lambda x: x[0])

    print("Sort done!")

    return pept_index
    

def identification(mzml_entry, pep_index, list_length):

    if "MSn spectrum" in mzml_entry:
        
        for precursor in mzml_entry["precursorList"]["precursor"]:

            if precursor["isolationWindow"]["ms level"] == 1:
                #auxiliary.print_tree(mzml_entry)
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
                        upper_index = list_length - 1
                    else:
                        upper_index = upper_index - 1

                    pep_index_slice = pep_index[lower_index:upper_index]

                    mzml_mz_array = mzml_entry["m/z array"]
                    mzml_intensity_array = mzml_entry["intensity array"]

                    binned_mzml_spectrum = binning(mzml_mz_array, mzml_intensity_array)

                    for pepts in pep_index_slice:

                        for pep in pepts[1]:

                            fasta_mz_array = np.array(sorted(list(fragments(pep, maxcharge=5)), key = float))
                            #fasta_mz_array = np.sort(np.array(list(fragments(pep, maxcharge=5)), key = float))
                    
                            binned_fasta_spectrum = binning(fasta_mz_array, theo_spect=True)

                            cross_corr = np.correlate(binned_mzml_spectrum, binned_fasta_spectrum, mode="full")

                            max_val = cross_corr.max()
                            second_max_val = np.partition(cross_corr.flatten(), -2)[-2]

                            xcorr_score = max_val - second_max_val

                            return xcorr_score
                        
    return None



def main(sample_filename : str, protein_database : str, processes : int, spectra_amount : int):
    start = time.time()
    
    with Path(protein_database).open("r", encoding="utf-8") as fasta_content:
        pept_list = create_pept_index(fasta_content)
    
    print("Computing Scores...")

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




