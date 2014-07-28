__author__ = 'dtgillis'
import bisect
import os


class Window():

    def __init__(self):

        self.window_start = None
        self.window_end = None
        self.window_depth_data = []

    def window_length(self):
        return self.window_end - self.window_start


class Interval():

    def __init__(self, chrm, interval_start, interval_end):

        self.chrm = chrm
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.read_count = []
        self.depth_of_coverage = []

    def __lt__(self, other):
        return self.interval_start < other.interval_start


class IntervalList():

    def __init__(self, interval_file):
        self.interval_list = []
        self.interval_file = interval_file

    def build_interval_list(self):

        for line in open(self.interval_file, 'r'):
            line = line.strip(os.linesep)
            if len(line) > 0:
                # construct interval
                chrm = line.split(':')[0]
                interval_start = int(line.split(':')[1].split('-')[0])
                interval_end = int(line.split(':')[1].split('-')[1])
                interval = Interval(chrm, interval_start, interval_end)
                # add interval to list
                self.interval_list.append(interval)
        # sort the list as a sanity check
        self.interval_list.sort()

    def add_interval(self, interval):
        """
        If a interval is added manually there will be a sort after so that the list
        stays in sorted order for fast interval access.
        :param interval:
        :return: None
        """
        self.interval_list.append(interval)
        self.interval_list.sort()

    def get_lower_bound(self, lower_bp):
        """
        This function is used to find the target region that
        contains this lower bound or is directly above it if it
        falls in a non target region
        Uses binary search of a sorted list
        :param lower_bp:
        :return: index of the interval that contains this lower bound or is
        above this lower bound if it does not fall in a target region
        """
        tmp_interval = Interval()
        tmp_interval.interval_start = lower_bp
        tmp_lower = bisect.bisect_left(self.interval_list, tmp_interval)
        #check to see if we have some overlap with the interval below
        #this happens if the interval below is like
        # interval.interval_start < lower_bp but interval.interval_end > lower_bp
        # interval [--------|------] [-----] where | is lower_bp
        interval_check = self.interval_list[tmp_lower - 1]
        if interval_check.interval_end >= lower_bp:
            return tmp_lower - 1
        else:
            return tmp_lower

    def get_upper_bound(self, upper_bp):
        """
        This function will return the upperbound location in the list
        If used for slicing you will need to add 1 to the loaction of
        the upper interval if it is to be included.
        Uses binary search
        :param upper_bp:
        :return: interval location
        """
        tmp_interval = Interval()
        tmp_interval.interval_start = upper_bp
        return bisect.bisect_left(self.interval_list, tmp_interval) - 1

    def number_of_intervals(self):
        return len(self.interval_list)


class CNVCall():
    """
    Class is representitive cnv for a region and sample
    """
    def __init__(self, window_start, window_end, cnv_state, sample):

        self.window_start = window_start
        self.window_end = window_end
        self.cnv_state = cnv_state
        self.sample = sample


from cnvpy.samtools_utils.io import SamFileDictionary


class DepthData():
    def __init__(self, bam_file_list, interval_list, binary=True):

        self.bam_files = SamFileDictionary(bam_file_list, binary=binary, open_files=True)
        self.intervals = IntervalList(interval_list)
        self.intervals.build_interval_list()























