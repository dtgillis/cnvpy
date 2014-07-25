__author__ = 'dtgillis'

import sys
import getopt
from cnvpy.depth_coverage.depth_parsers import SamToolsDepthParser
import cnvpy.depth_coverage.math_utils as doc_math
import cnvpy.depth_coverage.io_utils as io_utils
import os


def print_help():
    print """ here comes some help"""


def depth_using_sam_tools():
    interval_list = None
    out_dir = None
    bam_file_list = None
    window_size = None
    write_interval_stats = 0

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(
            argv, "hi:b:o:w:", ["interval_list=", "bam_file_list=", "out_dir=", "window_size="])
    except getopt.GetoptError:
        print 'cnv_caller.py -b <bamfilelist> -i <interval_list> -o <outputdir>'
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            print_help()
            sys.exit()
        elif opt in ("-i", "--interval_list"):
            interval_list = arg
        elif opt in ("-o", "--out_dir"):
            out_dir = arg
        elif opt in ("-b", "--bam_file_list"):
            bam_file_list = arg
        elif opt in ("-w", "--window_size"):
            window_size = int(arg)
        elif opt in ("--write_interval_stats",):
            write_interval_stats = arg
        else:
            assert False, "invalid option"

    assert interval_list is not None , "Specify interval list -i <interval_list>"
    assert out_dir is not None, "Specifiy output directory -o <out_dir>"
    assert bam_file_list is not None, "Specify bam file list -b <bam_file_list>"
    assert window_size is not None, "Specify a starting window size -w <window_size>"

    main_out_file_prefix = out_dir
    if not os.path.exists(main_out_file_prefix):
        try:
            os.makedirs(main_out_file_prefix)
        except OSError:
            print "unable to make output dir " + main_out_file_prefix
            sys.exit(1)

    depth_parser = SamToolsDepthParser(bam_file_list, interval_list, out_dir, window_size)

    for chrm in depth_parser.depth_data.interval_list.interval_dictionary:

        sys.stdout.write("Starting cnv_caller for " + chrm + os.linesep)
        depth_parser.get_interval_stats_by_sample(chrm, write_stats=(write_interval_stats == '1'))
        depth_parser.get_windows_in_intervals_by_chrm(chrm)
        io_utils.write_out_window_stats_by_chrm(depth_parser, chrm, main_out_file_prefix)
        depth_parser.clear_interval_information_for_chrm(chrm)

if __name__ == '__main__':

    depth_using_sam_tools()


