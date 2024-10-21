from pyteomics import mass
import math
from numpy import correlate

def masstocharge_to_dalton(mz : float, charge : int ):
    """ 
    Converts mass to charge to dalton
    """
    return mz * charge - 1.00794 * charge


def tolerance_bounds(ref_mass : float):
    """ 
    Returns the lower and upper tolerance bounds of the reference mass,
    as a tuple
    """

    lower_bound = ref_mass - ref_mass * 0.00001
    upper_bound = ref_mass + ref_mass * 0.00001
    return (lower_bound, upper_bound)


def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(
                            peptide[i:], ion_type=ion_type, charge=charge)
                    

def binary_search(pep_index, searched_mass, list_length):
    """ 
    Returns the index of the first mass in the list that is >= searched mass, returns -1 if 
    all list entrys are < searched mass
    """

    start = 0
    end = list_length - 1
    mid = 0

    if searched_mass > pep_index[end][0]:
        return -1
    
    while start <= end :

        mid = (start + end) // 2

        if pep_index[mid][0] >= searched_mass and pep_index[mid - 1][0] < searched_mass:
            break
        
        if pep_index[mid][0] < searched_mass:
            start = mid + 1

        elif pep_index[mid][0] > searched_mass:
            end = mid - 1
        
    return mid 


test_list = [[0.14], [0.4598], [0.78888], [1.9807], [2.333333], [6.65948503], [58.5432345], [543.3454324543], 
            [600.9999], [980.00009], [999.999], [1020.009], [1234.9089], [8999.8777], [9000.999]]
             

#print(binary_search(test_list, 0.009, len(test_list)))






