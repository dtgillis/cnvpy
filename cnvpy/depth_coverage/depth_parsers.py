__author__ = 'dtgillis'
import os
import sys
import pysam
from cnvpy.util_objects.depth_objects import DepthData


class BamFileListParser():

    def __init__(self, bam_list_file, binary=True):

        self.bam_list_file = bam_list_file
        self.bam_file_list = [line.strip().strip(os.sep) for line in open(self.bam_list_file, 'r')]
        self.binary = binary


class SamToolsDepthParser():

    def __init__(self, bam_list_file, interval_file, output_dir):
        self.depth_data = DepthData(bam_list_file, interval_file)
        self.output_dir = output_dir

    def parse_bam_files(self):
        sample_list = self.depth_data.bam_files.sample_list
        num_intervals = self.depth_data.intervals.number_of_intervals()
        sample_num = 0
        for sample in sample_list:
            sys.stdout.write("Parsing bam file for sample {0:s}\n".format(sample))
            sys.stdout.write("Parsing {0:d} Intervals\n".format(num_intervals))
            interval_num = 0
            samfp = self.depth_data.bam_files.sam_file_dict[sample]
            for interval in self.depth_data.intervals.interval_list:
                if interval_num % 1000 == 0:
                    sys.stdout.write("\n")
                elif interval_num % 100 == 0:
                    sys.stdout.write(".")
                self.get_interval_depth_of_coverage(samfp, interval, sample_num)
                interval_num += 1
            sample_num += 1
        self.depth_data.intervals.calculate_overall_avg()

        sys.stdout.write("\n Finished getting interval level stats\n\n")
        sys.stdout.write("\n------------------------------------------\n")

    def get_interval_depth_of_coverage(self, bam_file, interval, sample_num):
        # set up stats dictionary
        """
        Get all the depth stuff for the interval and then append it to the
        depth data of interval
        :param bam_file:
        :param interval:
        :param sample_name:
        """
        reads = bam_file.fetch(interval.chrm, interval.interval_start, interval.interval_end)
        read_count = 0
        for read in reads:

            assert isinstance(read, pysam.AlignedRead)
            if not read.is_duplicate and read.is_proper_pair:
                read_count += 1
        interval.read_count[sample_num] = read_count
        bases = bam_file.pileup(interval.chrm, interval.interval_start, interval.interval_end)
        for base in bases:

            assert isinstance(base, pysam.PileupProxy)
            if interval.interval_start <= base.pos <= interval.interval_end:

                interval.depth_of_coverage[sample_num, base.pos-interval.interval_start] = base.n
