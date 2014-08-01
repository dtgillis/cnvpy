from __future__ import division

__author__ = 'dtgillis'
import bisect
import os
from cnvpy.samtools_utils.io import SamFileDictionary
import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt


class Window():

    def __init__(self, window_start, window_end, depth_data, sample_num, cn_state):

        self.window_start = window_start
        self.window_end = window_end
        self.window_depth_data = depth_data
        self.num_target_regions = 1
        self.sample_num = sample_num
        self.window_done = False
        self.cn_state = cn_state

    def window_length(self):
        return self.window_end - self.window_start

    def add_interval(self, interval):

        if interval.interval_start < self.window_start:
            self.window_start = interval.interval_start
            self.window_depth_data = np.hstack((interval.depth_of_coverage, self.window_depth_data))
        else:
            self.window_end = interval.interval_end
            self.window_depth_data = np.hstack((self.window_depth_data, interval.depth_of_coverage))
        self.num_target_regions += 1

    def make_final_cnv_call(self, chrm_means):

        sample_data = self.window_depth_data[self.sample_num]
        probs = np.zeros_like(chrm_means)
        window_means = self.window_depth_data.mean(axis=1)
        target_constant = window_means.sum() / chrm_means.sum()

        for i in range(6):
            efficiency = np.floor(.5 * target_constant * chrm_means)
            probs[i] = np.nan_to_num(poisson.pmf(efficiency[self.sample_num] * i, sample_data.mean()))

        if probs.sum() != 0:
            overall_probs = probs / probs.sum()
            cnv_call = np.nonzero(probs.max() == probs)[0][0]
        else:
            cnv_call = 2

        return CNVCall(self.window_start, self.window_end, cnv_call, self.sample_num, self.num_target_regions)


class Interval():

    def __init__(self, chrm, interval_start, interval_end, num_samples):

        self.chrm = chrm
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.read_count = np.zeros(num_samples)
        if self.interval_start is None:
            self.depth_of_coverage = None
        else:
            self.depth_of_coverage = np.zeros([num_samples, interval_end - interval_start + 1])

    def __lt__(self, other):
        return self.interval_start < other.interval_start

    def get_interval_length(self):
        return self.interval_end - self.interval_start + 1

    def get_interval_avg(self):
        return self.depth_of_coverage.mean(axis=1)


class IntervalList():

    def __init__(self, interval_file, number_of_samples):
        self.interval_list = []
        self.interval_file = interval_file
        self.num_samples = number_of_samples
        self.overall_average = None

    def build_interval_list(self):

        for line in open(self.interval_file, 'r'):
            line = line.strip(os.linesep)
            if len(line) > 0:
                # construct interval
                chrm = line.split(':')[0]
                interval_start = int(line.split(':')[1].split('-')[0])
                interval_end = int(line.split(':')[1].split('-')[1])
                interval = Interval(chrm, interval_start, interval_end, self.num_samples)
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
        tmp_interval = Interval(None, None, None, None)
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
        tmp_interval = Interval(None, None, None, None)
        tmp_interval.interval_start = upper_bp
        return bisect.bisect_left(self.interval_list, tmp_interval)

    def number_of_intervals(self):
        return len(self.interval_list)

    def calculate_overall_avg(self):
        avg = np.zeros(self.num_samples)
        for interval in self.interval_list:
            avg += interval.get_interval_avg()
        self.overall_average = avg/self.number_of_intervals()


class CNVCall():
    """
    Class is representitive cnv for a region and sample
    """
    def __init__(self, window_start, window_end, cnv_state, sample, num_regions):

        self.window_start = window_start
        self.window_end = window_end
        self.cnv_state = cnv_state
        self.sample = sample
        self.num_regions = num_regions


class DepthData():
    """
    Main class Depth Data holds all information about the samples and
    intervals.

    contains IntervalList object and SamFileDictionary object.

    """
    def __init__(self, bam_file_list, interval_list, binary=True):

        self.bam_files = SamFileDictionary(bam_file_list, binary=binary, open_files=True)
        self.intervals = IntervalList(interval_list, len(self.bam_files.sample_list))
        self.intervals.build_interval_list()
        self.pass_one = []
        self.window_regions = []
        self.cnv_calls = []

    def target_calculations(self, mode='poisson', plotqq=False):

        # Calculate probs with target region as windows
        chrm_means = self.intervals.overall_average
        probs = []
        nb_probs = []
           # mass = importr('MASS')
        if len(self.pass_one) == 0:
            for interval in self.intervals.interval_list:
                window_means = interval.depth_of_coverage.mean(axis=1)
                target_constant = window_means.sum() / chrm_means.sum()
                efficiency = np.floor(target_constant * chrm_means)
                out_probs = np.zeros_like(efficiency)
                for i in range(self.intervals.num_samples):
                    #this prob / highest prob
                    if window_means[i] != 0:
                        out_probs[i] = \
                            poisson.pmf(int(efficiency[i]), window_means[i]) /\
                            poisson.pmf(int(window_means[i]), window_means[i])
                    else:
                        #TODO account for this with zero inflated model?
                        out_probs[i] = 0.0
                    if out_probs[i] < .05/22 and window_means.mean() > 8.0:
                        cn_state = window_means[i]/efficiency[i] * 2.0
                        self.window_regions.append(
                            Window(interval.interval_start,
                                   interval.interval_end,
                                   interval.depth_of_coverage,
                                   i, cn_state))

                probs.append(np.nan_to_num(out_probs))

            self.pass_one = np.array(probs)
            if plotqq:
                for i in range(self.pass_one.shape[1]):

                    y = -np.log(self.pass_one[:, i])
                    y = np.sort(y)
                    x = np.sort(-np.log(np.linspace(0, 1, self.pass_one.shape[0])))
                    plt.plot(x, y, '.')
                    plt.plot(x, x, '-')
                    plt.show()

    def window_calculations(self):

        windows_done = 0
        chrm_means = self.intervals.overall_average
        while windows_done != len(self.window_regions):

            # extend windows
            for window in self.window_regions:
                if not window.window_done:
                    loc = self.intervals.get_upper_bound(window.window_end)
                    if loc < self.intervals.number_of_intervals():
                        tmp_interval = self.intervals.interval_list[loc]
                        window.add_interval(tmp_interval)
                    else:
                        window.window_done = True
                        windows_done += 1

            for window in self.window_regions:

                if not window.window_done:
                    window_means = window.window_depth_data.mean(axis=1)
                    target_constant = window_means.sum() / chrm_means.sum()
                    efficiency = np.floor(window.cn_state * .5 * target_constant * chrm_means)
                    out_probs = np.zeros_like(efficiency)
                    sample = window.sample_num
                    #this prob / highest prob
                    if window_means[sample] != 0:
                        out_probs[sample] = \
                            poisson.pmf(int(efficiency[sample]), window_means[sample]) /\
                            poisson.pmf(int(window_means[sample]), window_means[sample])
                    else:
                        #TODO account for this with zero inflated model?
                        out_probs[sample] = 0.0
                    if out_probs[sample] < .95:
                        window.window_done = True
                        windows_done += 1
                        self.cnv_calls.append(window.make_final_cnv_call(chrm_means))