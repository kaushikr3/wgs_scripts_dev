import unittest

from python_analysis_scripts.generate_variant_outfiles import generate_annotation_series


class VCF_formatting(unittest.TestCase):
    # test vcf formatting and thresholding(?) here
    def test_something(self):
        self.assertEqual(True, False)


class GenBank_parsing_for_annotation(unittest.TestCase):
    gbk_path = ""
    gb_rec = [rec for rec in SeqIO.parse(gbk_path, "genbank")][0]
    gb_feats = [feat for feat in gb_rec.features if ((feat.type != "gene") & (feat.type != "source"))]

    start_coord_dict = {i: int(feat.location.start) for i, feat in enumerate(gb_feats)}
    end_coord_dict = {i: int(feat.location.end) for i, feat in enumerate(gb_feats)}

    def intergenic_feature_annotation(self):
        intergenic_annotation = generate_annotation_series()
        expected_output = []
        self.assertEqual(intergenic_annotation.tolist(), expected_output)

    def intergenic_end_of_genome_feature_annotation(self):

    def intergenic_beginning_of_genome_feature_annotation(self):

    def variant_in_final_coordinate_of_feature_annotation(self):

    def gene_feature_annotation(self):

    def CDS_feature_annotation(self):

    def nonCDS_feature_annotation(self):

    def multiple_CDS_feature_annotation(self):

    def CDS_and_nonCDS_feature_annotation(self):

    def frameshift_insertion_feature_annotation(self):

    def in_frame_insertion_feature_annotation(self):

    def frameshift_deletion_feature_annotation(self):

    def in_frame_deletion_feature_annotation(self):

    def silent_mutation_feature_annotation(self):

    def missense_mutation_feature_annotation(self):

    def truncate_mutation_feature_annotation(self):

class Annotation_formatting(unittest.TestCase):
    ## ensure that data merges into vcf data properly?

class MyTestCase(unittest.TestCase):
    # test lab variant cross-referencing here
    def test_something(self):
        self.assertEqual(True, False)




if __name__ == '__main__':
    unittest.main()
