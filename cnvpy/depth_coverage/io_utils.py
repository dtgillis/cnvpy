__author__ = 'dtgillis'

from cnvpy.depth_coverage.depth_parsers import SamToolsDepthParser
import os
from matplotlib.mlab import PCA
import pickle
import numpy as np
import cnvpy.depth_coverage.math_utils as doc_math
def write_out_window_stats_by_chrm(depth_parser, chrm, out_dir):

    assert isinstance(depth_parser, SamToolsDepthParser)
    interval_list = depth_parser.depth_data.interval_list.interval_dictionary[chrm]

    main_file = open(out_dir + os.sep + chrm + 'window_stats.dat', 'w')
    main_file.write("window_size")
    for sample in depth_parser.depth_data.bam_files.sam_file_dict.keys():
        main_file.write(",{0:s}_var,{0:s}_mean".format(sample))
    main_file.write(os.linesep)
    window_list = []
    average_chrm = doc_math.get_chrm_depth_avg_by_sample(interval_list)

    for interval in interval_list:
#        average_chrm = doc_math.get_chrm_depth_avg_by_sample(interval_list)
        for window in interval.windows:
            mid_window = []
            main_file.write("{0:d}".format(window.window_end - window.window_start + 1))
            for sample in average_chrm:
                main_file.write(",{0:f},{1:f}".format(window.window_var[sample], window.window_mean[sample]))
                mid_window.append(np.array(window.window_depth_by_sample[sample])[:, 1])

            main_file.write(os.linesep)
            window_list.append(mid_window)

    pickle.dump(window_list, open('/home/dtgillis/sim_capture/pickles/window_depth.p', 'wb'))
    avg = []
    for sample in average_chrm:
        avg.append(average_chrm[sample])

    pickle.dump(avg, open('/home/dtgillis/sim_capture/pickles/avg.p', 'wb'))
