# std import
from collections import defaultdict
import multiprocessing
from pathlib import Path
import time
from typing import DefaultDict, Set, List, TextIO, Tuple
import re
import matplotlib.pyplot as plt
import datetime

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

    Parameters
    ----------
    fasta_comtent : TextIO
        Opened fasta file
    
    Returns
    -------
    list
        Sorted eptide index list
    """

    print("Creating Peptide Index...")
    
    pept_index: DefaultDict[float, Set[str]] = defaultdict(set)

    variable_mods = {"[oxid]" : ["M"]}
    fixed_mods = {"[carb]" : ["C"]}

    mass.std_aa_comp['[oxid]'] = mass.Composition({'O': 1})
    mass.std_aa_comp['[carb]'] = mass.Composition({'H': 3, 'C' : 2, 'N' : 1, 'O' : 1})

    with read_fasta(fasta_content) as fasta:
        for (_, sequence) in fasta:
            peptides = parser.cleave(sequence, "trypsin", missed_cleavages=MAX_MISSED_CLEAVAGES, min_length=PEPTIDE_MIN_LENGTH, max_length=PEPTIDE_MAX_LENGTH)
            for peptide in peptides:
                if "X" not in peptide:

                    for modpept in parser.isoforms(peptide, variable_mods=variable_mods, fixed_mods=fixed_mods, max_mods=3):
                        pept_index[mass.calculate_mass(modpept)].add(modpept)


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
    Calculates scores and returns the identified result

    Parameters
    ----------
    mzml_entry : MzML spectrum
        Generated by the MzML reader iterator

    pep_index : list
        The peptide index list to search for

    list_length : int
        Length of the peptide index

    predict_spect : bool
        If spectrum prediction is used, parameter is set to true
    
    Returns
    -------
    list | None
        Returns the calculated results, if nothing was found, none is returned.
        
    """

    if mzml_entry["ms level"] == 2:

        spect_id = mzml_entry["id"]
        scan = int(spect_id.split("scan=")[1])

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

                if predict_spect:
                    binned_mzml_spectrum = binning(mzml_mz_array, mzml_intensity_array, spect_type=1)
                else:
                    binned_mzml_spectrum = binning(mzml_mz_array, mzml_intensity_array, spect_type=2)

                mzml_bins_size = binned_mzml_spectrum.size 
                
                #Append shift size on one of the arrays for shifted correlation
                shifted_binned_mzml_spectrum = np.insert(binned_mzml_spectrum, 0, np.zeros(int(SHIFT / 0.02)))
                shifted_binned_mzml_spectrum = np.append(shifted_binned_mzml_spectrum, np.zeros(int(SHIFT / 0.02)))

                xcorr_scores = []
                
                for pepts in pep_index_slice:

                    calc_neutral_mass = np.round(pepts[0], 6)

                    for pep in pepts[1]:

                        total_ions = 0

                        if predict_spect:

                            fasta_mz_array, fasta_int_array = predict_spectrum(pep, charge)
                            total_ions = fasta_mz_array.size
                            binned_fasta_spectrum = binning(fasta_mz_array, fasta_int_array, spect_type=1)
                            
                        else:
                            
                            fasta_mz_array = np.array(sorted(list(fragments(pep, maxcharge=min(charge-1, 3))), key = float))
                            total_ions = fasta_mz_array.size
                            binned_fasta_spectrum = binning(fasta_mz_array)

                        fasta_bins_size = binned_fasta_spectrum.size

                        if mzml_bins_size < fasta_bins_size:

                            binned_fasta_spectrum = binned_fasta_spectrum[:mzml_bins_size]

                        elif fasta_bins_size < mzml_bins_size:

                            binned_fasta_spectrum = np.append(binned_fasta_spectrum, np.zeros(mzml_bins_size-fasta_bins_size))

                        corr = np.correlate(shifted_binned_mzml_spectrum, binned_fasta_spectrum, "valid")
                        
                        # if scan in [71120]:

                        #     plt.figure(dpi=1200)
                        #     plt.plot(corr, linewidth=0.07, color='b')    
                        #     plt.title(f'Scan {scan} corr array ')    
                        #     plt.xlabel("Index") 
                        #     plt.ylabel("Correlation")                  
                        #     plt.savefig(f'Plots/corr_{scan}_ps={predict_spect}.png')


                        zeroshift_corr = corr[(corr.size // 2)] #Similarity at 0 offset
                        corr = np.delete(corr, (corr.size // 2)) #Delete similarity on Shift=0 before calculating background similarity
                        mean_corr = np.mean(corr) #Background similarity
                        
                        if predict_spect:
                            xcorr_score = np.round((zeroshift_corr - mean_corr) / 1000, 4) #Xcorr score
                        else:
                            xcorr_score = np.round((zeroshift_corr - mean_corr) / 10000, 4) 

                        matches = 0

                        for mzmlbin, fastabin in zip(binned_mzml_spectrum, binned_fasta_spectrum): #Count matches

                            if mzmlbin > 0 and fastabin > 0:
                                matches += 1

                        result = [scan, charge, np.round(mass_mzml, 6), calc_neutral_mass, xcorr_score, matches, total_ions,  pep] 

                        # if scan in [130051, 129688, 71120, 138707]:

                        #     plt.figure(dpi=1200)
                        #     plt.plot(binned_mzml_spectrum, linewidth=0.03, color='b')    
                        #     plt.plot(np.negative(binned_fasta_spectrum), linewidth=0.03, color='r')
                        #     plt.title(f'Scan: {scan} Score: {xcorr_score}')
                        #     plt.xlabel("Binned m/z")
                        #     plt.ylabel("Intensity")                       
                        #     plt.savefig(f'Plots/scan_{scan}_ps={predict_spect}.png')

                        xcorr_scores.append(result)
                    
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

    smallindex = pd.read_table("smallindex.txt", sep=' ')
    scanlist = [scan for scan in smallindex['scan']]

    with multiprocessing.Pool(processes) as pool, open(f'Results/{sample_filename.split('.')[0]}_ps={predict_spect}_result.tsv', "w") as outfile:

        print("scan\tcharge\texp_neutral_mass\tcalc_neutral_mass\txcorr\tions_matched\tions_total\tpeptide", file=outfile)
        
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
                        print(*r, file=outfile, sep='\t')

    end = time.time()
    print("Execution Time: ", str(datetime.timedelta(seconds= round(end - start))))




