__author__ = 'dtgillis'

import numpy as np

from cnvpy.util_objects.depth_objects import Window


class WindowStacker():
    """
    this works on a chrmosome basis
    """
    def __init__(self, depth_dict, window_size, chrm_interval_list, chrm):

        self.chrm_depth_dict = depth_dict
        self.window_size = window_size
        self.window_list = []
        self.chrm_interval_list = chrm_interval_list
        self.chrm = chrm


    def build_windows_from_chrm_dict(self):

        """
        :return:
        """
        bp_array, depth_data = self.chrm_depth_dict.create_numpy_objects_from_depth_data(self.chrm)

        for interval in self.chrm_interval_list:

            interval_length = interval.interval_end - interval.interval_start

            full_intervals = interval_length / self.window_size

            tmp_start = interval.interval_start

            for i in range(full_intervals):

                tmp_start += (i * self.window_size)

                tmp_index = np.nonzero(np.greater_equal(bp_array, tmp_start))

                start_index = tmp_index[0][0]

                tmp_start = bp_array[start_index]

                tmp_end = tmp_start + self.window_size - 1

                tmp_index = np.nonzero(np.less_equal(bp_array, tmp_end))

                if len(tmp_index) == 0:
                    end_index = bp_array[-1]
                else:
                    end_index = tmp_index[0][-1]
                tmp_window = Window()

                tmp_window.window_bp = bp_array[start_index:end_index+1]
                tmp_window.window_depth_data = depth_data[start_index:end_index+1]

                tmp_window.window_start = tmp_window.window_bp[0]
                tmp_window.window_end = tmp_window.window_bp[-1]

                self.window_list.append(tmp_window)
