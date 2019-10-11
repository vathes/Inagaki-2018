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
import tqdm

from . import reference, utilities, acquisition, analysis
from . import extracellular_path

schema = dj.schema(dj.config['custom'].get('database.prefix', '') + 'extracellular')
sess_data_dir = extracellular_path


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
        return NotImplementedError


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
    unit_depth: float  # (um)
    spike_waveform: longblob  # waveform(s) of each spike at each spike time (spike_time x waveform_timestamps)
    """

    def make(self, key):
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)

        if sess_data_file is None:
            raise FileNotFoundError(f'Extracellular import failed: ({key["subject_id"]} - {key["session_time"]})')

        mat_units = sio.loadmat(sess_data_file, struct_as_record=False, squeeze_me=True)['unit']
        for unit_idx, unit in tqdm.tqdm(enumerate(mat_units)):
            unit_key = dict(key,
                            unit_id=unit_idx,
                            channel_id=unit.channel,
                            unit_spike_width=unit.SpikeWidth,
                            unit_depth=unit.Depth,
                            spike_times=unit.SpikeTimes,
                            spike_waveform=unit.Spike_shpe_info.SpikeShape)
            self.insert1(unit_key, allow_direct_insert = True)


@schema
class TrialSegmentedUnitSpikeTimes(dj.Imported):
    definition = """
    -> UnitSpikeTimes
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_spike_times: longblob  # (s) with respect to the start of the trial
    """

    key_source = ProbeInsertion * analysis.TrialSegmentationSetting

    def make(self, key):
        # get data
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)

        if sess_data_file is None:
            raise FileNotFoundError(f'Extracellular import failed: ({key["subject_id"]} - {key["session_time"]})')

        mat_units = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['unit']

        unit_ids, spike_times = (UnitSpikeTimes & key).fetch('unit_id', 'spike_times')  # spike_times from all units
        trial_idx_of_spike = [mat_units[unit_id].Trial_idx_of_spike for unit_id in unit_ids]

        trial_keys = (acquisition.TrialSet.Trial & key).fetch('KEY')

        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')

        for trial_key in tqdm.tqdm(trial_keys):
            # get event time - this is done per trial
            try:
                event_time_point = analysis.get_event_time(event_name, trial_key)
            except analysis.EventChoiceError as e:
                print(f'Trial segmentation error - Msg: {str(e)}', file=sys.stderr)
                continue

            pre_stim_dur = float(pre_stim_dur)
            post_stim_dur = float(post_stim_dur)
            # check if pre/post stim dur is within start/stop time
            trial_start, trial_stop = (acquisition.TrialSet.Trial & trial_key).fetch1('start_time', 'stop_time')
            if trial_start and event_time_point - pre_stim_dur < 0:
                print('Warning: Out of bound prestimulus duration, select spikes from start-time (t=0)')
                pre_stim_dur = event_time_point
            if trial_stop and event_time_point + post_stim_dur > trial_stop:
                print('Warning: Out of bound poststimulus duration, set to trial end time')
                post_stim_dur = trial_stop - event_time_point

            # get spike times of each unit for this trial
            trial_spike_times = [spk[tr_spk_idx == trial_key['trial_id']]
                                 for spk, tr_spk_idx in zip(spike_times, trial_idx_of_spike)]
            # further segment based on pre/post stimulus duration
            seg_trial_spike_times = [spk[np.logical_and((spk >= (event_time_point - pre_stim_dur)),
                                                        (spk <= (event_time_point + post_stim_dur)))] - event_time_point
                                     for spk in trial_spike_times]

            self.insert(dict({**key, **trial_key},
                             unit_id=u_id,
                             segmented_spike_times=spk)
                        for u_id, spk in zip(unit_ids, seg_trial_spike_times))
