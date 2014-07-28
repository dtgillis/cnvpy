__author__ = 'dtgillis'
import sys

import numpy as np
from scipy.stats import poisson
from cnvpy.util_objects.depth_objects import CNVCall as cnv


def simple_bays_calc(sample_prob_matrix_dict, intervals):
    #TODO make this output some type of composite object instead of dictionary
    out_probs = dict()

    for sample in sample_prob_matrix_dict:

        #TODO this should be a try catch type deal
        p_c = sample_prob_matrix_dict[sample].mean(axis=0)

        p_c = np.nan_to_num(p_c)

        p_y = np.dot(p_c, sample_prob_matrix_dict[sample].T)

        p_y = np.nan_to_num(p_y)

        num = np.multiply(p_c, sample_prob_matrix_dict[sample])

        num = np.nan_to_num(num)

        out_probs[sample] = np.divide(num.T, p_y).T

        for i in range(out_probs[sample].shape[0]):

            out_probs[sample][i] /= out_probs[sample][i].sum()

    out_probs['coords'] = []

    for interval in intervals:

        for window in interval.windows:

            out_probs['coords'].append((window.window_start, window.window_end))

    return out_probs


def get_chrm_depth_avg_by_sample(chrm_interval_list):

    average_depth_by_sample = dict()

    for interval in chrm_interval_list:
        for window in interval.windows:
            window_depth = []
            for sample in sorted(window.window_depth_by_sample.keys()):
                if sample not in average_depth_by_sample:
                    average_depth_by_sample[sample] = .0
                tmp_depth = np.array(window.window_depth_by_sample[sample])
                window_depth.append(tmp_depth[:, 1])
                average_depth_by_sample[sample] += np.mean(tmp_depth[:, 1])
            window.window_depth_data = window_depth
            window.window_bp = tmp_depth[:, 0]

    return average_depth_by_sample


def calculate_window_stats(chrm_interval_list, chrm_means, average_depth_by_sample):
    sys.stdout.write("Calculating window level stats for " + str(len(chrm_interval_list)) + " intervals")
    copy_number = np.array([10**-8, 1, 2, 3, 4, 5, 6])
    interval_num = 0
    for interval in chrm_interval_list:
        if interval_num % 100 == 0:
            sys.stdout.write(".")
        if interval_num % 1000 == 0:
            sys.stdout.write("\n")
        interval_means = []
        for sample in interval.stats:
            interval_stats = interval.stats[sample]
            interval_rc = interval_stats['total_reads'] * 100 / float((interval.interval_end - interval.interval_start))
            interval_means.append(interval_rc)

        interval_means = np.array(interval_means)
        for window in interval.windows:
            window_means = np.mean(window.window_depth_data, axis=1)
            for sample in average_depth_by_sample:
                #interval_depth = interval.stats[sample]['total_reads'] * 100 / float(
                #    (interval.interval_end-interval.interval_start))
                target_eff_consant = \
                    (average_depth_by_sample[sample] * np.sum(window_means)) / (2.0 * np.sum(chrm_means))
                tmp_depth = int(np.array(window.window_depth_by_sample[sample])[:, 1].mean())
                tar_efficiency = copy_number * target_eff_consant

                for c in range(7):
                    if sample not in window.stats:
                        window.stats[sample] = []
                    frozen_pois = poisson(tar_efficiency[c])
                    window.stats[sample].append(frozen_pois.pmf(tmp_depth))

            for sample in window.stats:
                #sum all pvalues
                window.stats[sample] = np.array(window.stats[sample])
                window.stats[sample] /= window.stats[sample].sum()
        interval_num += 1

    sys.stdout.write("\nFinished windows stats.\n\n")
    sys.stdout.write("---------------------------------------\n")


def make_cnv_monte_carlo_call(prob_obj):

    cnv_state_dict = dict()
    for sample in prob_obj:

        if sample is not 'coords':

            cnv_state_dict[sample] = []

            for w_idx in range(len(prob_obj['coords'])):

                rand = np.random.random()

                cnv_vec = prob_obj[sample][w_idx]

                prb = .0
                cnv_state = 0
                while rand > prb:

                    prb += cnv_vec[0, cnv_state]
                    cnv_state += 1

                cnv_state_dict[sample].append(
                    cnv(prob_obj['coords'][w_idx][0], prob_obj['coords'][w_idx][1], cnv_state -1))

    return cnv_state_dict


def simple_cnv_call(prob_obj):

    cnv_state_dict = dict()
    for sample in prob_obj:

        if sample is not 'coords':

            cnv_state_dict[sample] = []

            for w_idx in range(len(prob_obj['coords'])):

                cnv_vec = prob_obj[sample][w_idx, :]
                if cnv_vec.max() > .999 and np.nonzero(cnv_vec == cnv_vec.max())[0] != 2:
                    cnv_state = np.nonzero(cnv_vec == cnv_vec.max())[0]
                    cnv_state_dict[sample].append(
                        cnv(prob_obj['coords'][w_idx][0], prob_obj['coords'][w_idx][1], cnv_state))

    return cnv_state_dict


def calculate_poisson_window_prob_for_chrm(chrm_interval_list):

    #hold onto sample depth data
    average_depth_by_sample = get_chrm_depth_avg_by_sample(chrm_interval_list)
    chrm_means = []
    for sample in sorted(average_depth_by_sample.keys()):
        average_depth_by_sample[sample] /= float(len(chrm_interval_list))
        chrm_means.append(average_depth_by_sample[sample])
    calculate_window_stats(chrm_interval_list, chrm_means, average_depth_by_sample)

    #TODO: make this a little more efficient
    sample_prob_matrices = dict()

    for interval in chrm_interval_list:
        for window in interval.windows:
            for sample in window.window_depth_by_sample:
                if sample not in sample_prob_matrices:
                    sample_prob_matrices[sample] = []
                sample_prob_matrices[sample].append(window.stats[sample])

    for sample in sample_prob_matrices:
        sample_prob_matrices[sample] = np.nan_to_num(np.array(sample_prob_matrices[sample]))

    out_probs = simple_bays_calc(sample_prob_matrices, chrm_interval_list)

    return simple_cnv_call(out_probs)
