import numpy as np
import math


def binning(mz_array, intensity_array=None, theo_spect=False, bin_width=0.02):
    """
    Bins spectra
    """

    #bins_filled = np.zeros(math.ceil(mz_array[mz_array.size - 1] / bin_width) + 1)
    bins_filled = np.zeros(int(2000/bin_width) + 1)

    if theo_spect:

        for mass in mz_array:
            if 200 <= mass <= 2000:

                index = int(mass // bin_width)
                bins_filled[index] = 1

                if index - 1 != -1:
                    bins_filled[index - 1] = max(bins_filled[index - 1], 0.5)

                bins_filled[index + 1] = max(bins_filled[index + 1], 0.5)
        

    else:
        
        for idx, mass in enumerate(mz_array):
            if  mass < 200 or mass > 2000:
                intensity_array[idx] = 0.0

        #Normalised intensity array 0-1
        intensity_array = intensity_array  / np.max(intensity_array)

        #Get top 100 intensities
        top_hundred_intensities = -np.sort(-intensity_array)[:min(100, intensity_array.size)]

        for mass, intensity in zip(mz_array, intensity_array):
            
            #Check if mass between 200 and 2000 and belongs to the top 100 intensities
            if intensity in top_hundred_intensities:

                index = int(mass // bin_width)
                bins_filled[index] = max(bins_filled[index], intensity)
        
    
    return bins_filled

#mz_array = np.array([56, 199, 200.01, 250, 300, 400, 700, 2000])
# int_array = np.array([2, 3, 6, 8, 9, 2, 90, 4, 56, 35])

#print(binning(mz_array, theo_spect=True)[10000])

