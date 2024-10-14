from pyteomics import mass

def masstocharge_to_dalton(mz : float, charge : int ):
    return mz * charge - 1.00794 * charge


def tolerance_bounds(ref_mass : float):
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

