__author__ = 'dtgillis'

import sys
import getopt
from cnvpy.depth_coverage.depth_parsers import SamToolsDepthParser
import cnvpy.depth_coverage.math_utils as doc_math
import os
import numpy as np

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
            argv, "hi:b:o:w:", ["interval_list=", "bam_file_list=", "out_dir=", "window_size=", "write_interval_stats="])
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

    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except OSError:
            print "error making path {0}".format(out_dir)

    main_out_file_prefix = out_dir + os.sep + 'cnv_calls.'

    depth_parser = SamToolsDepthParser(bam_file_list, interval_list, out_dir, window_size)

    sys.stdout.write("Starting cnv_caller for {0:d} samples".format(len(depth_parser.depth_data.bam_files.sample_list)))

    depth_parser.parse_bam_files()

    depth_parser.depth_data.target_calculations()

    depth_parser.depth_data.window_calculations()

    out_file_pointer = open(main_out_file_prefix + 'chrm1.chrm.out', 'w')
    #header stuff
    chrm = depth_parser.depth_data.intervals.interval_list[0].chrm
    sys.stdout.write("Writing cnv calls for " + chrm + " to file " + out_file_pointer.name + os.linesep)
    out_file_pointer.write('sample_name,chrm,num_exons,start,end,state,location' + os.linesep)
    sample_list = depth_parser.depth_data.bam_files.sample_list
    for cnv_call in depth_parser.depth_data.cnv_calls:
        if cnv_call.cnv_state != 2:
            out_file_pointer.write(
                "{0:s},{1:s},{5:d},{2:d},{3:d},{4:d},{1:s}:{2:d}-{3:d}\n".format(
                sample_list[cnv_call.sample], chrm,  int(cnv_call.window_start),
                int(cnv_call.window_end), int(cnv_call.cnv_state), cnv_call.num_regions))

    out_file_pointer.close()

#   depth_parser.clear_interval_information_for_chrm(chrm)

if __name__ == '__main__':

    depth_using_sam_tools()
