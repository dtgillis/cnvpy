__author__ = 'dtgillis'

import pysam
import os


class SamFileDictionary():

    def __init__(self, bam_files, binary=True, open_files=True):

        self.bam_files_list = bam_files
        self.binary = binary
        self.sam_file_dict = dict()
        self.sample_list = []

        for bam_file in open(self.bam_files_list, 'r'):

            bam_file = bam_file.strip(os.linesep)
            if self.binary:
                ext = ".bam"
                open_op = "rb"
            else:
                ext = ".sam"
                open_op = "r"

            sample_name = bam_file.split(os.sep)[-1].strip(ext)

            if sample_name in self.sam_file_dict:
                print "duplicate sample error"
            else:
                self.sam_file_dict[sample_name] = pysam.Samfile(bam_file, open_op)
                self.sample_list.append(sample_name)

    def close_bam_file(self, sample_name):

        if sample_name not in self.bam_files_list:
            print "error sample file not found"

        self.bam_files_list[sample_name].close()
















