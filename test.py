import unittest

from determine_hla import calculate_confidence
from seq2HLA import get_best_allele_per_group


class TestSeq2HLAMethods(unittest.TestCase):
    def test_get_best_allele_per_group(self):
        """
        Test that best alleles are selected per locus.
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
        #print(result)
        self.assertEqual(result, ({'A*01': 122, 'A*02': 0, 'B*45': 116, 'C*04': 74, 'G*29': 0}, {'A*01': 'A*01:01:02', 'B*45': 'B*45:07', 'C*04': 'C*04:02'}))

        
    def test_confidence_calculation_low_conf(self):
        """
        Test to make sure that confidence calculation is correct.
        """
        allele_max = 500
        allele_count_list = [500, 400, 300]
        result = calculate_confidence(allele_max, allele_count_list)
        print(result)
        self.assertEqual(result, 0.40444488206853557)


    def test_confidence_calculation_high_conf(self):
        """
        Test to make sure that confidence calculation is correct.
        """
        allele_max = 5000
        allele_count_list = [500, 400, 300]
        result = calculate_confidence(allele_max, allele_count_list)
        print(result)
        self.assertEqual(result, 0.0)

    def test_confidence_calculation_too_few_elements(self):
        """
        Test to make sure that too few elements in allele_count_list are not supported.
        """
        allele_max = 500
        allele_count_list = [500]
        result = calculate_confidence(allele_max, allele_count_list)
        self.assertEqual(result, 1.0)
        
if __name__ == '__main__':
    unittest.main()
