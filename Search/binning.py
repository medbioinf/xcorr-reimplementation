import numpy as np
import math


def binning(mz_array, intensity_array=None, spect_type=0, bin_width=0.02):
    """
    Bins spectra and returns their binned arrays.

    Parameters
    ----------
    mz_array : ndarray
        Array with m/z values

    intensity_array : ndarray
        Intensity array for the corresponding m/z array
    
    spect_type : int
        0 = Theoretical Spectrum,
        1 = Predicted Spectrum,
        2 = Experimental Spectrum

    bin_width : float
        Size of the bins for binning

    Returns
    -------
    ndarray
        The binned spectrum
    """

    if spect_type == 0:

        bins_filled = np.zeros(math.ceil(mz_array[mz_array.size - 1] / bin_width) + 1) 

        for mass in mz_array:

            index = int(mass // bin_width)
            bins_filled[index] = 50.0

            if index - 1 != -1:
                bins_filled[index - 1] = max(bins_filled[index - 1], 25.0)

            bins_filled[index + 1] = max(bins_filled[index + 1], 25.0)

        return bins_filled
    
    elif spect_type == 1:

        bins_filled = np.zeros(min(math.ceil(mz_array[mz_array.size - 1] / bin_width), int(2000 / 0.02)) + 1) 

        intensity_array = 50 * (intensity_array  / np.max(intensity_array))

        for mass, intensity in zip(mz_array, intensity_array):
            
            if mass <= 2000:
                
                index = int(mass // bin_width)
                bins_filled[index] = max(bins_filled[index], intensity)

        return bins_filled
    

    elif spect_type == 2:

        bins_filled = np.zeros(math.ceil(mz_array[mz_array.size - 1] / bin_width) + 1) 

        intensity_array = np.sqrt(intensity_array)

        for mass, intensity in zip(mz_array, intensity_array):
            
            index = int(mass // bin_width)
            bins_filled[index] = max(bins_filled[index], intensity)

        highest_ion = bins_filled.size 
        num_wins = 10
        win_size = int(highest_ion/num_wins) + 1

        norm_bins = np.array([]) 

        for i in range(0, len(bins_filled), win_size): 
            win = bins_filled[i:i + win_size]

            if np.max(win) != 0:
                win = 50 * (win  / np.max(win))

            norm_bins = np.append(norm_bins, win)
        
        
        del bins_filled
        #Get top 100 intensities
        # top_hundred_intensities = -np.sort(-norm_bins)[:min(100, norm_bins.size)]

        # for idx, intensity in enumerate(norm_bins):
        #     if intensity not in top_hundred_intensities:
        #         norm_bins[idx] = 0

        return norm_bins
    
    else:
        raise Exception("Not a valid spect_type")

