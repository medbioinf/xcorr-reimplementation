import numpy as np
import math

np.set_printoptions(threshold=np.inf, precision=10)

def binning(mz_array, intensity_array=None, theo_spect=False, bin_width=0.02):

    bins_filled = np.zeros(math.ceil(mz_array[len(mz_array) - 1] / bin_width) + 1)

    if theo_spect:

        for mass in mz_array:
            index = int(mass // bin_width)
            bins_filled[index] = 1

            if index - 1 != -1:
                bins_filled[index - 1] = max (bins_filled[index - 1], 0.5)

            bins_filled[index + 1] = max(bins_filled[index + 1], 0.5)

    else:
        #Normalised intensity array 0-1
        #intensity_array = (intensity_array - np.min(intensity_array)) / (np.max(intensity_array)-np.min(intensity_array)) 
        intensity_array = intensity_array  / np.max(intensity_array)
        
        for mass, intensity in zip(mz_array, intensity_array):

            index = int(mass // bin_width)
            bins_filled[index] = max(bins_filled[index], intensity)

    return bins_filled


