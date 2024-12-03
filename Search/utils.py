from pyteomics import mass, parser

def masstocharge_to_dalton(mz : float, charge : int ):
    """ 
    Converts mass to charge to dalton
    
    Parameters
    ----------
    mz : float
        Mass to charge value
    charge : int
        Charge state

    Returns
    -------
    float
        Returns mass in dalton
    """
    return (mz * charge) - (1.00794 * charge)


def tolerance_bounds(ref_mass : float):
    """ 
    Returns the lower and upper tolerance bounds of the reference mass,
    as a tuple

    Parameters
    ----------
    ref_mass : float
        The reference mass oround which the tolerance is calculated
    
    Returns
    -------
    (float, float)
        Tuple containing  lower and upper bound
    """
    #20 ppm
    lower_bound = ref_mass - ref_mass * 0.00002
    upper_bound = ref_mass + ref_mass * 0.00002
    return (lower_bound, upper_bound)


def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxcharge`.

    Parameters
    -----------
    peptide : str
        Pytemocs typed modX peptide string 

    types : tuple
        Types of ions for fragmentation.
        Default = ('b', 'y')

    maxcharge : int
        Maximum charge of the fragments

    Yields
    ------
    float
        Generates m/z
    """
    peptide = parser.parse(peptide)

    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.calculate_mass(peptide[:i], ion_type=ion_type, charge=charge)

                else:
                    yield mass.calculate_mass(peptide[i:], ion_type=ion_type, charge=charge)
                        

def binary_search(pep_index: list[tuple[float, tuple[str]]], searched_mass, list_length):
    """ 
    Returns the index of the first mass in the list that is >= searched mass, returns -1 if 
    all list entrys are < searched mass

    Parameters
    ----------
    pep_index : list
        The peptide index list to search for

    searched_mass : float
        The mass to search for
    
    list_length : int
        Length of the pep_index list

    Returns
    -------
    int
        Returns the index of the found entry in the pep_index.
        Returns -1 if nothing is found
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
