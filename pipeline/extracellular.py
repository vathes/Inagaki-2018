'''
Schema of extracellular information.
'''
import re
import os
import sys
from datetime import datetime

import numpy as np
import scipy.io as sio
import datajoint as dj
import h5py as h5
import tqdm

from . import reference, utilities, acquisition, analysis

schema = dj.schema(dj.config.get('database.prefix', '') + 'extracellular')


@schema
class ProbeInsertion(dj.Manual):
    definition = """ # Description of probe insertion details during extracellular recording
    -> acquisition.Session
    -> reference.Probe
    -> reference.BrainLocation
    """


@schema
class Voltage(dj.Imported):
    definition = """
    -> ProbeInsertion
    ---
    voltage: longblob   # (mV)
    voltage_start_time: float # (second) first timepoint of voltage recording
    voltage_sampling_rate: float # (Hz) sampling rate of voltage recording
    """

    def make(self, key):
        # this function implements the ingestion of raw extracellular data into the pipeline
        return None


@schema
class UnitSpikeTimes(dj.Imported):
    definition = """ 
    -> ProbeInsertion
    unit_id : smallint
    ---
    -> reference.Probe.Channel
    spike_times: longblob  # (s) time of each spike, with respect to the start of session 
    unit_cell_type='N/A': varchar(32)  # e.g. cell-type of this unit (e.g. wide width, narrow width spiking)
    unit_spike_width: float  # (ms) spike width of this unit, from bottom peak to next positive peak or time point spike terminates
    unit_depth: float  # (mm)
    spike_waveform: longblob  # waveform(s) of each spike at each spike time (spike_time x waveform_timestamps)
    """

    def make(self, key):
        sess_data_dir = os.path.join('.', 'data', 'SiliconProbeData')
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)

        if sess_data_file is None:
            print(f'Extracellular import failed: ({key["subject_id"]} - {key["session_time"]})', file=sys.stderr)
            return

        mat_units = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['unit']
        for unit_idx, unit in tqdm.tqdm(enumerate(mat_units)):
            unit_key = dict(key,
                            unit_id = unit_idx,
                            channel_id = unit.channel,
                            unit_spike_width = unit.SpikeWidth,
                            unit_depth = unit.Depth,
                            spike_times = unit.SpikeTimes,
                            spike_waveform = unit.Spike_shpe_info.SpikeShape)
            self.insert1(unit_key, allow_direct_insert = True)


@schema
class TrialSegmentedUnitSpikeTimes(dj.Imported):
    definition = """
    -> UnitSpikeTimes
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_spike_times: longblob
    """

    def make(self, key):
        # get data
        sess_data_dir = os.path.join('.', 'data', 'SiliconProbeData')
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)

        if sess_data_file is None:
            print(f'Extracellular import failed: ({key["subject_id"]} - {key["session_time"]})', file=sys.stderr)
            return
        mat_units = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['unit']

        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get event time
        try:
            event_time_point = analysis.get_event_time(event_name, key)
        except analysis.EventChoiceError as e:
            print(f'Trial segmentation error - Msg: {str(e)}')
            return

        pre_stim_dur = float(pre_stim_dur)
        post_stim_dur = float(post_stim_dur)
        # check if pre/post stim dur is within start/stop time
        trial_start, trial_stop = (acquisition.TrialSet.Trial & key).fetch1('start_time', 'stop_time')
        if trial_start and event_time_point - pre_stim_dur < trial_start:
            print('Warning: Out of bound prestimulus duration, set to 0')
            pre_stim_dur = 0
        if trial_stop and event_time_point + post_stim_dur > trial_stop:
            print('Warning: Out of bound poststimulus duration, set to trial end time')
            post_stim_dur = trial_stop - event_time_point

        # get raw & segment
        spike_times = (UnitSpikeTimes & key).fetch1('spike_times')[
            mat_units[key['unit_id']].Trial_idx_of_spike - 1 == key['trial_id']]  # -1 to account for MATLAB 1-based indexing
        key['segmented_spike_times'] = spike_times[np.logical_and(
            (spike_times >= (event_time_point - pre_stim_dur)),
            (spike_times <= (event_time_point + post_stim_dur)))] - event_time_point

        self.insert1(key)
        print(f'Perform trial-segmentation of spike times for unit: {key["unit_id"]} and trial: {key["trial_id"]}')
