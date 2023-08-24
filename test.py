import unittest

from seq2HLA import get_best_allele_per_group


class TestSeq2HLAMethods(unittest.TestCase):
    def test_get_best_allele_per_group(self):
        """
        Test that it can sum a list of integers
        """
        readspergroup = {
            'A*01': 0,
            'A*02': 0,
            'B*45': 0,
            'C*04': 0,
            'G*29': 0
        }
        readcount = {
            'A*01:05': 90,
            'A*01:01:02': 122,
            'B*45:07': 116,
            'C*04:02': 74
        }
        
        result = get_best_allele_per_group(readspergroup, readcount)
        print(result)
        self.assertEqual(result, ({'A*01': 122, 'A*02': 0, 'B*45': 116, 'C*04': 74, 'G*29': 0}, {'A*01': 'A*01:01:02', 'B*45': 'B*45:07', 'C*04': 'C*04:02'}))

        
if __name__ == '__main__':
    unittest.main()
