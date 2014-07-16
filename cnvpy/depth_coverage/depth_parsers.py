__author__ = 'dtgillis'
import os
import pickle
import sys
import math

import pysam


class IntervalParser():

    def __init__(self, interval_file):

        self.interval_file = interval_file
        self.interval_dictionary = dict()

    def build_interval_dictionary(self, window_size):

        from cnvpy.util_objects.depth_objects import Interval

        for line in open(self.interval_file, 'r'):

            interval = Interval(window_size)

            interval.chrm = line.split(':')[0]
            interval.interval_start = int(line.split(':')[1].split('-')[0])
            interval.interval_end = int(line.split(':')[1].split('-')[1].strip(os.sep))

            if interval.chrm not in self.interval_dictionary:

                self.interval_dictionary[interval.chrm] = []

            self.interval_dictionary[interval.chrm].append(interval)


#TODO: make the bam list file two columns (file) (binary)
# So we can mix the types
class BamFileListParser():

    def __init__(self, bam_list_file, binary=True):

        self.bam_list_file = bam_list_file
        self.bam_file_list = [line.strip().strip(os.sep) for line in open(self.bam_list_file, 'r')]
        self.binary = binary

from cnvpy.util_objects.depth_objects import DepthData
from cnvpy.util_objects.depth_objects import Window


class SamToolsDepthParser():

    def __init__(self, bam_list_file, interval_file, output_dir, window_size):
        self.depth_data = DepthData(bam_list_file, interval_file, window_size)
        self.output_dir = output_dir
        # add the stats dict for each sample we are looking at
        for chrm in self.depth_data.interval_list.interval_dictionary:
            for sample in self.depth_data.bam_files.sam_file_dict:
                chrm_interval_list = self.depth_data.interval_list.interval_dictionary[chrm]
                for interval in chrm_interval_list:
                    interval.stats[sample] = dict()

    def get_interval_stats_by_sample(self, chrm, write_stats=False):

        sample_list = self.depth_data.bam_files.sam_file_dict.keys()

        chrm_interval_list = self.depth_data.interval_list.interval_dictionary[chrm]
        interval_num = 0
        sys.stdout.write(
            "Creating interval stats for %s intervals" % len(chrm_interval_list))
        for interval in chrm_interval_list:
            if interval_num % 100 == 0:
                sys.stdout.write(".")
            if interval_num % 1000 == 0:
                sys.stdout.write("\n")

            for sample in sample_list:

                samfp = self.depth_data.bam_files.sam_file_dict[sample]

                self.get_interval_read_stats(samfp, interval, sample)
            interval_num += 1

        sys.stdout.write("\n Finished getting interval level stats\n\n")
        sys.stdout.write("\n------------------------------------------\n")
        if write_stats:

            for sample in sample_list:
                stats_out = open(self.output_dir + os.sep + sample + '.' + chrm + '.interval_stats.dat', 'w')
                stats_out.write('Int_size,reads,avg_insert,unpaired_percent,gc_content' + os.linesep)
                for interval in chrm_interval_list:
                    stats = interval.stats[sample]
                    stats_out.write(str(interval.interval_end-interval.interval_start))
                    for stat in ['total_reads', 'avg_insert_size', 'unpaired_percent', 'gc_content']:
                        stats_out.write(",%s" % stats[stat])
                    stats_out.write(os.linesep)
                stats_out.close()

            stats_out = open(self.output_dir + os.sep + chrm + '.interval_read_count.dat', 'w')
            stats_out.write('Chrm,Int_start,Int_end,Size')
            for sample in sample_list:
                stats_out.write(',' +sample)
            stats_out.write(os.linesep)

            for interval in chrm_interval_list:
                stats_out.write(
                    "%s,%s,%s,%s" % (chrm, interval.interval_start, interval.interval_end,
                                  interval.interval_end-interval.interval_start))
                for sample in sample_list:
                    stats_out.write(",%s" % interval.stats[sample]['total_reads'])
                stats_out.write(os.linesep)

    def get_interval_read_stats(self, bam_file, interval, sample_name):
        # set up stats dictionary
        """

        :param bam_file:
        :param interval:
        :param sample_name:
        """
        interval.stats[sample_name]['total_insert_size'] = 0
        interval.stats[sample_name]['total_reads'] = 0
        interval.stats[sample_name]['unpaired_reads'] = 0
        interval.stats[sample_name]['total_quality'] = 0.0
        interval.stats[sample_name]['total_g_c'] = 0
        interval.stats[sample_name]['total_a_t'] = 0
        reads = bam_file.fetch(interval.chrm, interval.interval_start, interval.interval_end)

        for read in reads:

            assert isinstance(read, pysam.AlignedRead)
            if not read.is_duplicate:

                if read.is_read1:

                    if read.is_proper_pair:
                        insert_size = math.fabs(read.tlen)
                        interval.stats[sample_name]['total_insert_size'] += insert_size
                    else:
                        interval.stats[sample_name]['unpaired_reads'] += 1
                seq_string = read.seq

                interval.stats[sample_name]['total_g_c'] += seq_string.count('G') + seq_string.count('C')
                interval.stats[sample_name]['total_a_t'] += seq_string.count('A') + seq_string.count('T')
                interval.stats[sample_name]['total_quality'] += read.mapq
                interval.stats[sample_name]['total_reads'] += 1

        if interval.stats[sample_name]['total_reads'] != 0:
            interval.stats[sample_name]['avg_insert_size'] = \
                interval.stats[sample_name]['total_insert_size'] / float(interval.stats[sample_name]['total_reads'])
            interval.stats[sample_name]['avg_quality'] = \
                interval.stats[sample_name]['total_quality'] / float(interval.stats[sample_name]['total_reads'])
            interval.stats[sample_name]['unpaired_percent'] = \
                interval.stats[sample_name]['unpaired_reads'] / float(interval.stats[sample_name]['total_reads'])
            interval.stats[sample_name]['gc_content'] = \
                interval.stats[sample_name]['total_g_c'] / float(interval.stats[sample_name]['total_g_c']
                                                                 + interval.stats[sample_name]['total_a_t'])
        else:
            interval.stats[sample_name]['avg_insert_size'] = 0.0
            interval.stats[sample_name]['avg_quality'] = 0.0
            interval.stats[sample_name]['unpaired_percent'] = 0.0
            interval.stats[sample_name]['gc_content'] = 0.0

    def get_windows_in_intervals_by_chrm(self, chrm):

        if chrm not in self.depth_data.interval_list.interval_dictionary:
            print "error in interval list"
        interval_num = 0
        sys.stdout.write(
            "Creating windows for %s intervals" % len(self.depth_data.interval_list.interval_dictionary[chrm]))
        for interval in self.depth_data.interval_list.interval_dictionary[chrm]:
            if interval_num % 100 == 0:
                sys.stdout.write(".")
            if interval_num % 1000 == 0:
                sys.stdout.write("\n")

            self.get_windows_in_interval_for_sample_list(
                    interval, interval.windows_size)
            interval_num += 1

        sys.stdout.write("\nFinished creating windows    \n\n")
        sys.stdout.write("-----------------------------------\n")

    def get_windows_in_interval_for_sample_list(self, interval, window_size):

        total_windows = (interval.interval_end - interval.interval_start) / window_size

        if total_windows == 0:
            total_windows = 1

        window_start = interval.interval_start

        for i in range(total_windows):

            tmp_window = Window()

            tmp_window.window_start = window_start

            if window_start + window_size < interval.interval_end and i != total_windows-1:
                window_end = window_start + window_size - 1
            else:
                window_end = interval.interval_end

            tmp_window.window_end = window_end

            for sample in self.depth_data.bam_files.sam_file_dict:

                bam_file = self.depth_data.bam_files.sam_file_dict[sample]
                assert isinstance(bam_file, pysam.Samfile)

                bases = bam_file.pileup(interval.chrm, tmp_window.window_start, tmp_window.window_end)

                tmp_depth_list = []
                # this is a pointer to basepair
                # should help track regions with
                # no reads
                tmp_cur_base = tmp_window.window_start

                for base in bases:

                    assert isinstance(base, pysam.PileupProxy)

                    if tmp_window.window_start <= base.pos <= tmp_window.window_end:

                        if tmp_cur_base != base.pos:
                            while tmp_cur_base < base.pos:
                                #pad with 0s beggining or in side
                                tmp_depth_list.append((tmp_cur_base, 0))
                                tmp_cur_base += 1

                        tmp_depth_list.append((base.pos, base.n))
                        tmp_cur_base += 1

                if tmp_cur_base <= tmp_window.window_end:
                    # stuff at the end with no read...
                    while tmp_cur_base <= tmp_window.window_end:
                        #pad with 0s
                        tmp_depth_list.append((tmp_cur_base, 0))
                        tmp_cur_base += 1

                tmp_window.window_depth_by_sample[sample] = tmp_depth_list

            interval.windows.append(tmp_window)
            window_start += window_size

    def pickle_intervals_for_chrm(self):
        pass
        # for chrm in self.depth_data.interval_list.interval_dictionary:
        #
        #     chrm_intervals = self.depth_data.interval_list.interval_dictionary[chrm]
        #     pickle.dump(chrm_intervals, open('' + chrm + '.interval_depth.p', 'wb'))

    def clear_interval_information_for_chrm(self, chrm):

        self.depth_data.interval_list.interval_dictionary[chrm] = None