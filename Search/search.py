# std import
from collections import defaultdict
import multiprocessing
from pathlib import Path
import time
from typing import DefaultDict, Set, List, TextIO, Tuple
import re
import matplotlib.pyplot as plt

# 3rd party import
from pyteomics import mzml, parser, auxiliary, mass
from pyteomics.fasta import read as read_fasta
import numpy as np
import pandas as pd

# internal imports
from Search.utils import fragments, masstocharge_to_dalton, tolerance_bounds, binary_search
from Search.binning import binning
from Search.predicted_spect import predict_spectrum

import pickle

PEPTIDE_MIN_LENGTH = 6
PEPTIDE_MAX_LENGTH = 50
MAX_MISSED_CLEAVAGES = 2
SHIFT = 75

def create_pept_index(fasta_content: TextIO) -> List[Tuple[float, Tuple[str]]]:
    """ 
    Creates the sorted peptide index List, contains the mass and peptides
    """

    print("Creating Peptide Index...")
    
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
    

def identification(mzml_entry, pep_index, list_length, predict_spect, scanlist):
    """ 
    Returns the identified result as [scan, xcorrscore, peptide]
    """

    if mzml_entry["ms level"] == 2:
        #auxiliary.print_tree(mzml_entry)
        spect_id = mzml_entry["id"]
        scan = int(spect_id.split("scan=")[1])
        #scan = int(re.search("scan=\d+", spect_id).group(0)[5:])

        #if scan == 71120:
        #if scan in [71120, 131926]:
        #if scan in scanlist:
        #if scan == 131926:
        if scan in [71120, 102968, 113974, 102326, 70466, 107413, 121180, 103621, 124639, 131926, 114319, 131364, 87329, 107790]:

            for precursor in mzml_entry["precursorList"]["precursor"]:

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

                    pep_index_slice = pep_index[lower_index:upper_index]

                    mzml_mz_array = mzml_entry["m/z array"]
                    mzml_intensity_array = mzml_entry["intensity array"]

                    binned_mzml_spectrum = binning(mzml_mz_array, mzml_intensity_array)

                    #Get mzml spect before appending shiftsize for correct plotting
                    binned_mzml_spectrum_before_shift = binned_mzml_spectrum #numpy.copy()?

                    mzml_bins_size = binned_mzml_spectrum.size 
                    
                    #Append shift size on one of the arrays for shift
                    binned_mzml_spectrum = np.insert(binned_mzml_spectrum, 0, np.zeros(int(SHIFT / 0.02)))
                    binned_mzml_spectrum = np.append(binned_mzml_spectrum, np.zeros(int(SHIFT / 0.02)))

                    xcorr_scores = []

                    for pepts in pep_index_slice:

                        for pep in pepts[1]:

                            if predict_spect:

                                fasta_mz_array, fasta_int_array = predict_spectrum(pep, charge-1)
                                binned_fasta_spectrum = binning(fasta_mz_array, fasta_int_array)

                            else:
                                fasta_mz_array = np.array(sorted(list(fragments(pep, maxcharge=charge-1)), key = float))
                                binned_fasta_spectrum = binning(fasta_mz_array, theo_spect=True)
    
                            fasta_bins_size = binned_fasta_spectrum.size

                            if mzml_bins_size < fasta_bins_size:

                                binned_fasta_spectrum = binned_fasta_spectrum[:mzml_bins_size]

                            elif fasta_bins_size < mzml_bins_size:

                                binned_fasta_spectrum = np.append(binned_fasta_spectrum, np.zeros(mzml_bins_size-fasta_bins_size))
            

                            corr = np.correlate(binned_mzml_spectrum, binned_fasta_spectrum, "valid")

                            zeroshift_corr = corr[(corr.size // 2)] #Similarity at 0 offset
                            corr = np.delete(corr, (corr.size // 2)) #Delete similarity on Shift=0 before calculating background similarity
                            mean_corr = np.mean(corr) #Background similarity
                            
                            xcorr_score = zeroshift_corr - mean_corr #Xcorr score

                            matches = 0

                            for mzmlbin, fastabin in zip(binned_mzml_spectrum_before_shift, binned_fasta_spectrum): #Count matches

                                if mzmlbin > 0 and fastabin > 0:
                                    matches += 1

                            result = [scan, xcorr_score, matches, pep] 
                            print(result)
                            xcorr_scores.append(result)

                            # with open(f"testdata/binned_mzml_spectrum_{scan}.pkl", "bw") as f:
                            #     pickle.dump(binned_mzml_spectrum_before_shift, f)

                            #if scan in [100440, 131364, 131926, 132489, 138755]:

                            # plt.figure(dpi=1200)
                            # plt.plot(binned_mzml_spectrum_before_shift, linewidth=0.03, color='b')    
                            # plt.plot(np.negative(binned_fasta_spectrum), linewidth=0.03, color='r')
                            # plt.title(f'Scan: {scan} Score: {xcorr_score}')                       

                            # if predict_spect:

                            #     plt.savefig(f'Plots/scan_{scan}_ps.png')

                            # else:
                                
                            #     plt.savefig(f'Plots/scan_{scan}_comet_top.png')

                    return xcorr_scores
                        
    return None



def main(sample_filename : str, protein_database : str, processes : int, spectra_amount : int, predict_spect : bool):
    start = time.time()
    
    with Path(protein_database).open("r", encoding="utf-8") as fasta_content:
        pept_list = create_pept_index(fasta_content)
    
    print("Computing Scores...")

    list_length = len(pept_list)

    manager = multiprocessing.Manager()
    pep_index = manager.list(pept_list) # <- peptide index

    total_results = []

    smallindex = pd.read_table("smallindex.txt", sep=' ')
    scanlist = [scan for scan in smallindex['scan']]

    with multiprocessing.Pool(processes) as pool, open(f'{sample_filename.split('.')[0]}_ps={predict_spect}_result.txt', "w") as outfile:

        print("Run against: ", protein_database, file=outfile)
        print("Scan XcorrScore Matches Peptide", file=outfile)
        
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
                pool.apply_async(identification, args=(spec, pep_index, list_length, predict_spect, scanlist))
                for spec in spec_buffer
            ]

            for result in results:
                res = result.get()

                if res is not None:
                    
                    for r in res:
                        print(*r, file=outfile)

    end = time.time()
    print("Execution Time: ", end - start)




