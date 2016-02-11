__author__ = 'kupa'

import unittest
import os
import subprocess


class CarnivoresTestPhyRe(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Used to run PhyRe on the sample list running everything else"""

        self.sample_file = 'tests/carnivores_sample.txt'
        self.sample_out = 'tests/carnivores_sample.out'
        self.funnel_file = 'tests/carnivores_funnel.out'
        self.master_file = 'tests/carnivores_list.txt'
        self.d1 = 40
        self.d2 = 100
        self.missing_data = 'y'
        self.permutations = '100'

        subprocess.call('python PhyRe.py '+self.sample_file+' '+self.master_file+' '+str(self.d1)+' '+str(self.d2)+
                        ' -m '+self.missing_data+' -p '+self.permutations, shell=True)


    def test_funnel_line_number(self):
        funnel_output_length = 0
        with open(self.funnel_file, 'r') as fp:
            funnel_output_length = len(fp.readlines())
        self.assertEqual(funnel_output_length, 67,
                         'Wrong Funnel File Length')


    def test_sample_out_line_number(self):
        sample_file_output_length = 0
        with open(self.sample_out, 'r') as fp:
            sample_file_output_length = len(fp.readlines())
        self.assertEqual(sample_file_output_length, 29,
                         'Wrong Sample File Length')

    def test_taxa_number_and_path_length (self):
        """Number of taxa and path lengths for each taxonomic level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        taxons = ['family', 'subfamily', 'genus', 'species']
        taxons_number = [11, 23, 129, 271]
        path_lengths = [32, 18, 29, 18]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and path lengths for each taxonomic level:') == 0:
                break
        i = i+1

        #parsing strings of the output file from 4th till 8th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertGreater(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] - 1, 'Error in '+taxons[tu]+ ' path length ')
            self.assertLess(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] + 1, 'Error in '+taxons[tu]+ ' path length ')

    def test_taxa_number_and_pairwise_comparison_number (self):
        """Number of taxa and pairwise comparisons  at each taxon level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        taxons = ['family', 'subfamily', 'genus', 'species']
        taxons_number = [11, 25, 60, 72]
        pairwise_comparisons = [4268, 490, 298, 56]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and pairwise comparisons') == 0:
                break
        i = i+1

        #parsing strings of the output file from 15th till 19th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertEqual(int(lines[i+tu].split('\t')[2]) , pairwise_comparisons[tu], 'Error in '+taxons[tu]+ ' pairwise comparison ')

    def test_key_metrics_block (self):
        """Testing Key Metrics"""

        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        metrics = ['Average taxonomic distinctness', 'Variation in taxonomic distinctness',
                  'Minimum taxonomic distinctness', 'Maximum taxonomic distinctness', 'von Eulers index of imbalance']
        values = [92, 280, 66, 96, 0]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find(metrics[0]) == 0:
                break

        #parsing strings of the output file from 22th till 26th
        for tu in range(len(metrics)):
            self.assertGreater(float(lines[i+tu].split('=')[1]), values[tu] - 1, 'Error in '+metrics[tu])
            self.assertLess(float(lines[i+tu].split('=')[1]), values[tu] + 1, 'Error in '+metrics[tu])


    @classmethod
    def tearDownClass(self):
        """called once, after all tests, if setUpClass successful"""

        os.remove(self.sample_out)
        os.remove(self.funnel_file)
        self.sample_file = None
        self.funnel_file = None
        self.sample_out = None
        self.master_file = None
        self.missing_data = None
        self.permutations = None


class BivalvesTestPhyRe(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Used to run PhyRe on the sample list running everything else"""

        self.sample_file = 'tests/bivalves_sample.txt'
        self.sample_out = 'tests/bivalves_sample.out'
        self.funnel_file = 'tests/bivalves_funnel.out'
        self.master_file = 'tests/bivalves_list.txt'
        self.d1 = 5
        self.d2 = 15
        self.missing_data = 'y'
        self.permutations = '100'

        subprocess.call('python PhyRe.py '+self.sample_file+' '+self.master_file+' '+str(self.d1)+' '+str(self.d2)+
                        ' -m '+self.missing_data+' -p '+self.permutations, shell=True)


    def test_funnel_line_number(self):
        funnel_output_length = 0
        with open(self.funnel_file, 'r') as fp:
            funnel_output_length = len(fp.readlines())
        self.assertEqual(funnel_output_length, 17,
                         'Wrong Funnel File Length')


    def test_sample_out_line_number(self):
        sample_file_output_length = 0
        with open(self.sample_out, 'r') as fp:
            sample_file_output_length = len(fp.readlines())
        self.assertEqual(sample_file_output_length, 29,
                         'Wrong Sample File Length')

    def test_taxa_number_and_path_length (self):
        """Number of taxa and path lengths for each taxonomic level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        taxons = ['Subclass', 'Order', 'Family', 'Genus']
        taxons_number = [5, 17, 217, 3404]
        path_lengths = [23, 20, 27, 27]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and path lengths for each taxonomic level:') == 0:
                break
        i = i+1

        #parsing strings of the output file from 4th till 8th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertGreater(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] - 1, 'Error in '+taxons[tu]+ ' path length ')
            self.assertLess(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] + 1, 'Error in '+taxons[tu]+ ' path length ')

    def test_taxa_number_and_pairwise_comparison_number (self):
        """Number of taxa and pairwise comparisons  at each taxon level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        taxons = ['Subclass', 'Order', 'Family', 'Genus']
        taxons_number = [3, 5, 8, 9]
        pairwise_comparisons = [52, 10, 8, 2]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and pairwise comparisons') == 0:
                break
        i = i+1

        #parsing strings of the output file from 15th till 19th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertEqual(int(lines[i+tu].split('\t')[2]) , pairwise_comparisons[tu], 'Error in '+taxons[tu]+ ' pairwise comparison ')

    def test_key_metrics_block (self):
        """Testing Key Metrics"""

        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        metrics = ['Average taxonomic distinctness', 'Variation in taxonomic distinctness',
                  'Minimum taxonomic distinctness', 'Maximum taxonomic distinctness', 'von Eulers index of imbalance']
        values = [89, 340, 79, 90, 0]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find(metrics[0]) == 0:
                break

        #parsing strings of the output file from 22th till 26th
        for tu in range(len(metrics)):
            self.assertGreater(float(lines[i+tu].split('=')[1]), values[tu] - 1, 'Error in '+metrics[tu])
            self.assertLess(float(lines[i+tu].split('=')[1]), values[tu] + 1, 'Error in '+metrics[tu])


    @classmethod
    def tearDownClass(self):
        """called once, after all tests, if setUpClass successful"""

        os.remove(self.sample_out)
        os.remove(self.funnel_file)
        self.sample_file = None
        self.funnel_file = None
        self.sample_out = None
        self.master_file = None
        self.missing_data = None
        self.permutations = None


class ColeoidsTestPhyRe(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Used to run PhyRe on the sample list running everything else"""

        self.sample_file = 'tests/coleoids_sample.txt'
        self.sample_out = 'tests/coleoids_sample.out'
        self.funnel_file = 'tests/coleoids_funnel.out'
        self.master_file = 'tests/coleoids_list.txt'
        self.d1 = 20
        self.d2 = 40
        self.missing_data = 'y'
        self.permutations = '100'

        subprocess.call('python PhyRe.py '+self.sample_file+' '+self.master_file+' '+str(self.d1)+' '+str(self.d2)+
                        ' -m '+self.missing_data+' -p '+self.permutations, shell=True)


    def test_funnel_line_number(self):
        funnel_output_length = 0
        with open(self.funnel_file, 'r') as fp:
            funnel_output_length = len(fp.readlines())
        self.assertEqual(funnel_output_length, 27,
                         'Wrong Funnel File Length')


    def test_sample_out_line_number(self):
        sample_file_output_length = 0
        with open(self.sample_out, 'r') as fp:
            sample_file_output_length = len(fp.readlines())
        self.assertEqual(sample_file_output_length, 33,
                         'Wrong Sample File Length')

    def test_taxa_number_and_path_length (self):
        """Number of taxa and path lengths for each taxonomic level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        taxons = ['Order', 'Suborder', 'Superfamily', 'Family', 'Subfamily', 'Genus']
        taxons_number = [7, 7, 8, 52, 19, 220]
        path_lengths = [22, 0, 3, 22, 26, 24]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and path lengths for each taxonomic level:') == 0:
                break
        i = i+1

        #parsing strings of the output file from 4th till 8th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertGreater(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] - 1, 'Error in '+taxons[tu]+ ' path length ')
            self.assertLess(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] + 1, 'Error in '+taxons[tu]+ ' path length ')

    def test_taxa_number_and_pairwise_comparison_number (self):
        """Number of taxa and pairwise comparisons  at each taxon level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()

        taxons = ['Order', 'Suborder', 'Superfamily', 'Family', 'Subfamily', 'Genus']
        taxons_number = [4, 7, 9, 17, 21, 30]
        pairwise_comparisons = [614, 100, 40, 60, 24, 32]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and pairwise comparisons') == 0:
                break
        i = i+1

        #parsing strings of the output file from 15th till 19th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertEqual(int(lines[i+tu].split('\t')[2]) , pairwise_comparisons[tu], 'Error in '+taxons[tu]+ ' pairwise comparison ')

    def test_key_metrics_block (self):
        """Testing Key Metrics"""

        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        metrics = ['Average taxonomic distinctness', 'Variation in taxonomic distinctness',
                  'Minimum taxonomic distinctness', 'Maximum taxonomic distinctness', 'von Eulers index of imbalance']
        values = [90, 315, 72, 92, 0]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find(metrics[0]) == 0:
                break

        #parsing strings of the output file from 22th till 26th
        for tu in range(len(metrics)):
            self.assertGreater(float(lines[i+tu].split('=')[1]), values[tu] - 1, 'Error in '+metrics[tu])
            self.assertLess(float(lines[i+tu].split('=')[1]), values[tu] + 1, 'Error in '+metrics[tu])


    @classmethod
    def tearDownClass(self):
        """called once, after all tests, if setUpClass successful"""

        os.remove(self.sample_out)
        os.remove(self.funnel_file)
        self.sample_file = None
        self.funnel_file = None
        self.sample_out = None
        self.master_file = None
        self.missing_data = None
        self.permutations = None


class TermitesTestPhyRe(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Used to run PhyRe on the sample list running everything else"""

        self.sample_file = 'tests/termites_sample.txt'
        self.sample_out = 'tests/termites_sample.out'
        self.funnel_file = 'tests/termites_funnel.out'
        self.master_file = 'tests/termites_list.txt'
        self.d1 = 20
        self.d2 = 60
        self.missing_data = 'y'
        self.permutations = '100'

        subprocess.call('python PhyRe.py '+self.sample_file+' '+self.master_file+' '+str(self.d1)+' '+str(self.d2)+
                        ' -m '+self.missing_data+' -p '+self.permutations, shell=True)


    def test_funnel_line_number(self):
        funnel_output_length = 0
        with open(self.funnel_file, 'r') as fp:
            funnel_output_length = len(fp.readlines())
        self.assertEqual(funnel_output_length, 47,
                         'Wrong Funnel File Length')


    def test_sample_out_line_number(self):
        sample_file_output_length = 0
        with open(self.sample_out, 'r') as fp:
            sample_file_output_length = len(fp.readlines())
        self.assertEqual(sample_file_output_length, 33,
                         'Wrong Sample File Length')

    def test_taxa_number_and_path_length (self):
        """Number of taxa and path lengths for each taxonomic level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()

        taxons = ['family', 'subfamily', 'alliance', 'genus', 'section', 'species']
        taxons_number = [7, 17, 2, 282, 11, 2760]
        path_lengths = [15, 10, 18, 18, 18, 18]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and path lengths for each taxonomic level:') == 0:
                break
        i = i+1

        #parsing strings of the output file from 4th till 8th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertGreater(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] - 1, 'Error in '+taxons[tu]+ ' path length ')
            self.assertLess(float(lines[i+tu].split('\t')[2]) , path_lengths[tu] + 1, 'Error in '+taxons[tu]+ ' path length ')

    def test_taxa_number_and_pairwise_comparison_number (self):
        """Number of taxa and pairwise comparisons  at each taxon level test"""
        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        taxons = ['family', 'subfamily', 'alliance', 'genus', 'section', 'species']
        taxons_number = [7, 17, 17, 40, 40, 40]
        pairwise_comparisons = [1214, 206, 0, 140, 0, 0]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find('Number of taxa and pairwise comparisons') == 0:
                break
        i = i+1

        #parsing strings of the output file from 15th till 19th
        for tu in range(len(taxons)):
            self.assertEqual(int(lines[i+tu].split('\t')[1]), taxons_number[tu], 'Error in '+taxons[tu]+ ' number ')
            self.assertEqual(int(lines[i+tu].split('\t')[2]) , pairwise_comparisons[tu], 'Error in '+taxons[tu]+ ' pairwise comparison ')

    def test_key_metrics_block (self):
        """Testing Key Metrics"""

        with open(self.sample_out, 'r') as fp:
            lines = fp.readlines()


        metrics = ['Average taxonomic distinctness', 'Variation in taxonomic distinctness',
                  'Minimum taxonomic distinctness', 'Maximum taxonomic distinctness', 'von Eulers index of imbalance']
        values = [93, 177, 78, 96, 0]

        for i in range(len(lines)):
            if lines[i].split('\t')[0].find(metrics[0]) == 0:
                break

        #parsing strings of the output file from 22th till 26th
        for tu in range(len(metrics)):
            self.assertGreater(float(lines[i+tu].split('=')[1]), values[tu] - 1, 'Error in '+metrics[tu])
            self.assertLess(float(lines[i+tu].split('=')[1]), values[tu] + 1, 'Error in '+metrics[tu])


    @classmethod
    def tearDownClass(self):
        """called once, after all tests, if setUpClass successful"""

        os.remove(self.sample_out)
        os.remove(self.funnel_file)
        self.sample_file = None
        self.funnel_file = None
        self.sample_out = None
        self.master_file = None
        self.missing_data = None
        self.permutations = None


if __name__ == '__main__':
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(CarnivoresTestPhyRe))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(BivalvesTestPhyRe))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(ColeoidsTestPhyRe))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TermitesTestPhyRe))
    unittest.TextTestRunner(verbosity=2).run(suites)
