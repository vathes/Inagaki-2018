'''
Schema of intracellular information.
'''
import re
import os
import sys
from datetime import datetime

import numpy as np
import scipy.io as sio
import datajoint as dj

from . import reference, utilities, acquisition, analysis
from . import intracellular_path

schema = dj.schema(dj.config.get('database.prefix', '') + 'intracellular')
sess_data_dir = os.path.join(intracellular_path, 'Data')


@schema
class Cell(dj.Manual):
    definition = """ # A cell undergone intracellular recording in this session
    -> acquisition.Session
    cell_id: varchar(36) # a string identifying the cell in which this intracellular recording is concerning
    ---
    cell_type: enum('Pyr', 'GABA', 'excitatory', 'inhibitory', 'N/A')
    cell_depth: float  # (um)
    -> reference.BrainLocation
    -> reference.WholeCellDevice
    """


@schema
class MembranePotential(dj.Imported):
    definition = """
    -> Cell
    ---
    membrane_potential: longblob  # (mV) membrane potential recording at this cell
    membrane_potential_wo_spike: longblob # (mV) membrane potential without spike data, derived from membrane potential recording    
    membrane_potential_start_time: float # (s) first timepoint of membrane potential recording
    membrane_potential_sampling_rate: float # (Hz) sampling rate of membrane potential recording
    """

    def make(self, key):
        # ============ Dataset ============
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            raise FileNotFoundError(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})')

        mat_data = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['wholeCell']

        #  ============= Now read the data and start ingesting =============
        print(f'Insert membrane potential data for: {key["cell_id"]}')
        # -- MembranePotential
        self.insert1(dict(
            key,
            membrane_potential=mat_data.recording_data.Vm,
            membrane_potential_wo_spike=mat_data.recording_data.Vm_wo_spike,
            membrane_potential_start_time=0,
            membrane_potential_sampling_rate=mat_data.recording_data.sample_rate))


@schema
class CurrentInjection(dj.Imported):
    definition = """
    -> Cell
    ---
    injected_current: float  # (pA) amplitude setting of the current injection routine
    current_injection: longblob
    current_injection_start_time: float  # first timepoint of current injection recording
    current_injection_sampling_rate: float  # (Hz) sampling rate of current injection recording
    """
    
    # -- CurrentInjection - only available for EPSP session
    key_source = Cell & (acquisition.Session.ExperimentType & 'experiment_type = "EPSP"')

    def make(self, key):
        # ============ Dataset ============
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            raise FileNotFoundError(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})')

        mat_data = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['wholeCell']

        #  ============= Now read the data and start ingesting =============
        print(f'Insert current injection data for: {key["cell_id"]}')
        self.insert1(dict(
            key,
            injected_current=mat_data.meta_data.injected_current,
            current_injection=mat_data.recording_data.Output_700B,
            current_injection_start_time=0,
            current_injection_sampling_rate=mat_data.recording_data.sample_rate))


@schema
class CellSpikeTimes(dj.Imported):
    definition = """
    -> Cell
    ---
    spike_times: longblob  # (s) time of each spike, with respect to the start of session 
    """

    def make(self, key):
        # ============ Dataset ============
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            raise FileNotFoundError(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})')

        mat_data = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['wholeCell']

        #  ============= Now read the data and start ingesting =============
        print(f'Insert spikes data for: {key["cell_id"]}')
        # -- Spike
        self.insert1(dict(
            key,
            spike_times=mat_data.recording_data.spike_peak_bin / mat_data.recording_data.sample_rate))


@schema
class TrialSegmentedMembranePotential(dj.Computed):
    definition = """
    -> MembranePotential
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_mp: longblob   
    segmented_mp_wo_spike: longblob
    """

    key_source = MembranePotential * acquisition.TrialSet * analysis.TrialSegmentationSetting

    def make(self, key):
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get raw
        fs, first_time_point, Vm_wo_spike, Vm_w_spike = (MembranePotential & key).fetch1(
            'membrane_potential_sampling_rate', 'membrane_potential_start_time', 'membrane_potential_wo_spike',
            'membrane_potential')

        # Limit to insert size of 15 per insert
        insert_size = 15
        trial_lists = utilities.split_list((acquisition.TrialSet.Trial & key).fetch('KEY'), insert_size)

        for b_idx, trials in enumerate(trial_lists):
            self.insert(dict({**key, **trial_key},
                             segmented_mp = perform_trial_segmentation(trial_key, event_name,
                                                                        pre_stim_dur, post_stim_dur,
                                                                        Vm_w_spike, fs, first_time_point),
                             segmented_mp_wo_spike = perform_trial_segmentation(trial_key, event_name,
                                                                                 pre_stim_dur, post_stim_dur,
                                                                                 Vm_wo_spike, fs, first_time_point))
                        for trial_key in trials)


@schema
class TrialSegmentedCurrentInjection(dj.Computed):
    definition = """
    -> CurrentInjection
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_current_injection: longblob
    """

    key_source = CurrentInjection * acquisition.TrialSet * analysis.TrialSegmentationSetting

    def make(self, key):
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get raw
        fs, first_time_point, current_injection = (CurrentInjection & key).fetch1(
            'current_injection_sampling_rate', 'current_injection_start_time', 'current_injection')
        self.insert(dict({**key, **trial_key},
                         segmented_current_injection = perform_trial_segmentation(trial_key, event_name,
                                                                                  pre_stim_dur, post_stim_dur,
                                                                                  current_injection, fs,
                                                                                  first_time_point))
                    for trial_key in (acquisition.TrialSet.Trial & key).fetch('KEY'))


@schema
class TrialSegmentedCellSpikeTimes(dj.Computed):
    definition = """
    -> CellSpikeTimes
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_spike_times: longblob
    """

    def make(self, key):
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get event time
        try:
            event_time_point = analysis.get_event_time(event_name, key)
        except analysis.EventChoiceError as e:
            raise e

        pre_stim_dur = float(pre_stim_dur)
        post_stim_dur = float(post_stim_dur)
        # check if pre/post stim dur is within start/stop time
        trial_start, trial_stop = (acquisition.TrialSet.Trial & key).fetch1('start_time', 'stop_time')
        event_time_point = event_time_point + trial_start

        if trial_start and event_time_point - pre_stim_dur < 0:
            print('Warning: Out of bound pre-stim dur, select spikes from start-time (t=0)')
            pre_stim_dur = event_time_point - trial_start
        if event_time_point + post_stim_dur > trial_stop:
            print('Warning: Out of bound post-stim dur, set to trial end time')
            post_stim_dur = trial_stop - event_time_point

        # get raw & segment
        spike_times = (CellSpikeTimes & key).fetch1('spike_times')
        key['segmented_spike_times'] = spike_times[np.logical_and(
            (spike_times >= (event_time_point - pre_stim_dur)),
            (spike_times <= (event_time_point + post_stim_dur)))] - event_time_point

        self.insert1(key)
        print(f'Perform trial-seg spike times for cell: {key["cell_id"]} - trial: {key["trial_id"]}')


def perform_trial_segmentation(trial_key, event_name, pre_stim_dur, post_stim_dur, data, fs, first_time_point):
    # get event time
    try:
        event_time_point = analysis.get_event_time(event_name, trial_key)
    except analysis.EventChoiceError as e:
        raise e
    #
    pre_stim_dur = float(pre_stim_dur)
    post_stim_dur = float(post_stim_dur)
    sample_total = int((post_stim_dur + pre_stim_dur)*fs) + 1
    # check if pre/post stim dur is within start/stop time, if not, pad with NaNs
    trial_start, trial_stop = (acquisition.TrialSet.Trial & trial_key).fetch1('start_time', 'stop_time')
    event_time_point = event_time_point + trial_start

    pre_stim_nan_count = 0
    post_stim_nan_count = 0
    if trial_start and event_time_point - pre_stim_dur < trial_start:
        pre_stim_nan_count = int(np.ceil((trial_start - (event_time_point - pre_stim_dur)) * fs))
        pre_stim_dur = event_time_point - trial_start
        print(f'Warning: Out-of-bound pre-stim dur, pad {pre_stim_nan_count} NaNs')
    if trial_stop and event_time_point + post_stim_dur > trial_stop:
        post_stim_nan_count = int(np.ceil((event_time_point + post_stim_dur - trial_stop) * fs))
        post_stim_dur = trial_stop - event_time_point
        print(f'Warning: Out-of-bound post-stim dur, pad {post_stim_nan_count} NaNs')

    event_sample_point = (event_time_point - first_time_point) * fs
    sample_points_to_extract = range((event_sample_point - pre_stim_dur * fs).astype(int),
                                     (event_sample_point + post_stim_dur * fs + 1).astype(int))
    segmented_data = data[sample_points_to_extract]
    # pad with NaNs
    segmented_data = np.hstack((np.full(pre_stim_nan_count, np.nan), segmented_data,
                                np.full(post_stim_nan_count, np.nan)))
    return segmented_data[:sample_total]


