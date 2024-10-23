import numpy

numpy.set_printoptions(threshold=numpy.inf)

def binning(mz_array, intensity_array=None, bin_width=0.02, maxvalue=10000):

    bin_edges = []
    edge = 0.0

    while edge <= (mz_array[len(mz_array) - 1]):
        bin_edges.append(edge)
        edge = round(edge + bin_width, 4)

    bins_filled = [0 for i in range(len(bin_edges))]

    if intensity_array == None:

        for mass in mz_array:
            bins_filled[numpy.digitize(mass, bin_edges) - 1] = maxvalue

    else:
        for mass in mz_array:

            edge_index = numpy.digitize(mass, bin_edges) - 1
            int_array_index = mz_array.index(mass)
            if bins_filled[edge_index] < intensity_array[int_array_index]:
                bins_filled[edge_index] = intensity_array[int_array_index]

    return bins_filled

    

mz_array = [0.001, 0.14, 0.4598, 0.78888, 1.9808, 1.9807, 2.333333, 10.0, 11.00001]
int_array = [27, 0.4, 598, 9, 888, 8, 807, 10.9, 10.900001]


#print(binning(mz_array))