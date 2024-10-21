import numpy

numpy.set_printoptions(threshold=numpy.inf)

def binning(spectrum, bin_width=1.0005):

    bin_edges = []
    edge = 0.0

    while edge < spectrum[len(spectrum) - 1 ] + bin_width:
        bin_edges.append(edge)
        edge = round(edge + bin_width, 5)

    bins_filled = []
    for mass in spectrum:
        bins_filled.append(numpy.digitize(mass, bin_edges))
    
    return bin_edges, bins_filled

test_list = [0.14, 0.4598, 0.78888, 1.9807, 2.333333, 6.65948503, 58.5432345, 543.3454324543, 
            600.9999, 980.00009, 999.999, 1020.009, 1234.9089, 8999.8777, 9000.999]


# print(binning(test_list))