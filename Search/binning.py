import numpy as np
import math


def binning(mz_array, intensity_array=None, theo_spect=False, bin_width=0.02):
    """
    Bins spectra
    """

    bins_filled = np.zeros(math.ceil(mz_array[mz_array.size - 1] / bin_width) + 1)

    if theo_spect:

        for mass in mz_array:
            index = int(mass // bin_width)
            bins_filled[index] = 1

            if index - 1 != -1:
                bins_filled[index - 1] = max(bins_filled[index - 1], 0.5)

            bins_filled[index + 1] = max(bins_filled[index + 1], 0.5)

    else:
        #Normalised intensity array 0-1
        intensity_array = intensity_array  / np.max(intensity_array)

        #Get top 100 intensities
        top_hundred_intensities = -np.sort(-intensity_array)[:min(100, intensity_array.size)]

        for mass, intensity in zip(mz_array, intensity_array):
            
            #Check if mass between 200 and 2000 and belongs to the top 100 intensities
            if 200 <= mass <= 2000 and intensity in top_hundred_intensities:

                index = int(mass // bin_width)
                bins_filled[index] = max(bins_filled[index], intensity)

    return bins_filled


