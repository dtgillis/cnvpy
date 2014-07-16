__author__ = 'dtgillis'
import numpy as np

class Window():

    def __init__(self):

        self.window_start = None
        self.window_end = None
        self.window_depth_data = None
        self.window_depth_by_sample = dict()
        self.window_bp = None
        self.stats = dict()


class Interval():

    def __init__(self, window_size):

        self.chrm = None
        self.interval_start = None
        self.interval_end = None
        self.windows_size = window_size
        self.windows = []
        self.stats = dict()
        self.done = False
        self.windows_prbs = []


class CNVCall():

    def __init__(self, window_start, window_end, cnv_state):

        self.window_start = window_start
        self.window_end = window_end
        self.cnv_state = cnv_state


class DepthDataGATK():

    def __init__(self):

        self.main_dict = dict()

    def chrm_set(self, chrm):

        return chrm in self.main_dict

    def add_chrm(self, chrm):

        self.main_dict[chrm] = dict()
        self.main_dict[chrm]['bp'] = []
        self.main_dict[chrm]['depth'] = []

    def get_chrm(self, chrm):

        if chrm in self.main_dict:

            return self.main_dict[chrm]

    def add_depth_data(self, chrm, bp, depth_data):

        if chrm not in self.main_dict:

            print "error " + chrm + " not in dictionary"

        self.main_dict[chrm]['bp'].append(bp)
        self.main_dict[chrm]['depth'].append(depth_data)

    def create_numpy_objects_from_depth_data(self, chrm):
        """

        :param chrm:
        :return: tuple of numpy arrays basepairs, depth_data

        basepairs is 1 by p , p is number of basepairs in chrm
        depthdata is p by n , n is number of samples
        """
        np_depth = np.array(self.main_dict[chrm]['depth'], dtype=np.int)
        np_bp = np.array(self.main_dict[chrm]['bp'], dtype=np.int)

        return np_bp, np_depth

from cnvpy.samtools_utils.io import SamFileDictionary
from cnvpy.depth_coverage.depth_parsers import IntervalParser

class DepthData():
    def __init__(self, bam_file_list, interval_list, window_size, binary=True):

        self.bam_files = SamFileDictionary(bam_file_list, binary=binary, open_files=True)
        self.interval_list = IntervalParser(interval_list)
        self.interval_list.build_interval_dictionary(window_size)























