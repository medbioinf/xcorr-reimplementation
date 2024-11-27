import sys
import unittest

sys.path.append("../")

from Search.utils import tolerance_bounds, masstocharge_to_dalton, fragments, binary_search
from Search.binning import binning

class TestUtils(unittest.TestCase):

    def test_tolerance_bounds(self):
        lower, upper = tolerance_bounds(90.9)
        self.assertEqual(lower, 90.898182)
        self.assertEqual(upper, 90.901818)
    
    def test_masstocharge_to_dalton(self):
        pass

    def test_fragments(self):
        pass

    def test_binary_search(self):
        pass

class TestBinning(unittest.TestCase):

    def test_binning(self):
        pass

if __name__ == "__main__":
    unittest.main()