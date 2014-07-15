__author__ = 'dtgillis'

import sys
import getopt
from cnvpy.depth_coverage.depth_parsers import SamToolsDepthParser
import cnvpy.depth_coverage.math_utils as doc_math
import os


def print_help():
    pass


def depth_using_sam_tools():
    interval_list = None
    out_dir = None
    bam_file_list = None

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hi:b:o:", ["interval_list=", "bam_file_list=", "out_dir="])
    except getopt.GetoptError:
        print 'cnv_caller.py -b <bamfilelist> -i <interval_list> -o <outputdir>'
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print_help()
        elif opt in ("-i", "--interval_list"):
            interval_list = arg
        elif opt in ("-o", "--out_dir"):
            out_dir = arg
        elif opt in ("-b", "--bam_file_list"):
            bam_file_list = arg
        else:
            assert False, "invalid option"

    assert interval_list is not None , "Specify interval list -i <interval_list>"
    assert out_dir is not None, "Specifiy output directory -o <out_dir>"
    assert bam_file_list is not None, "Specify bam file list -b <bam_file_list>"

    main_out_file_prefix = out_dir + os.sep + 'cnv_calls.'

    depth_parser = SamToolsDepthParser(bam_file_list, interval_list)
    depth_parser.get_interval_stats_by_sample()
    for chrm in depth_parser.depth_data.interval_list.interval_dictionary:

        depth_parser.get_windows_in_intervals_by_chrm(chrm)
        sample_prob = doc_math.calculate_poisson_window_prob_for_chrm(
            depth_parser.depth_data.interval_list.interval_dictionary[chrm])
        out_file_pointer = open(main_out_file_prefix  + chrm + '.out', 'w')

        #header stuff

        out_file_pointer.write('sample_name,chrm,start,end,state' + os.linesep)
        for sample in sample_prob:
            for cnv_call in sample_prob[sample]:
                out_file_pointer.write(
                    "%s,%s,%s,%s,%s\n" % (sample, chrm, cnv_call.window_start,
                                          cnv_call.window_end, cnv_call.cnv_state[0]))
        out_file_pointer.close()

        depth_parser.clear_interval_information_for_chrm(chrm)

if __name__ == '__main__':

    depth_using_sam_tools()

