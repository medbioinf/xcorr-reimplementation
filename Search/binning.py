import numpy as np
import math

np.set_printoptions(threshold=np.inf, precision=10)

def binning(mz_array, intensity_array=None, theo_spect=False, bin_width=0.02):

    if theo_spect == True:

        bins_filled = np.zeros(math.ceil(mz_array[len(mz_array) - 1] / bin_width) + 2)

        for mass in mz_array:
            index = int(mass // bin_width) + 1
            bins_filled[index] = 1

            if bins_filled[index - 1] == 0:
                bins_filled[index - 1] = 0.5

            if bins_filled[index + 1] == 0:
                bins_filled[index + 1] = 0.5

    else:
        bins_filled = np.zeros(math.ceil(mz_array[len(mz_array) - 1] / bin_width))

        #Normalised intensity array 0-1
        intensity_array = (intensity_array - np.min(intensity_array)) / (np.max(intensity_array)-np.min(intensity_array))

        for mass, intensity in zip(mz_array, intensity_array):

            index = int(mass // bin_width)
            if bins_filled[index] < intensity:
                bins_filled[index] = intensity

    return bins_filled
    


mz_array = np.array([0.001, 0.14, 0.4598, 0.78888, 1.9808, 1.9807, 2.333333, 10.0])
# int_array = np.array([27, 0.4, 598, 9, 888, 8, 807, 10.9])

# # # # spect1 = binning(mz_array, int_array)
# # # # spect2 = binning(mz_array, theo_spect=1)

# # # # print(np.correlate(spect1, spect2, mode="full"))
# print(binning(mz_array, theo_spect=True))
# # print(np.round((mz_array - np.min(mz_array)) / (np.max(mz_array)-np.min(mz_array)), 8))




