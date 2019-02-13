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
import h5py as h5

from . import reference, utilities, acquisition, analysis

schema = dj.schema(dj.config.get('database.prefix', '') + 'intracellular')


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
    membrane_potential: longblob    # (mV) membrane potential recording at this cell
    membrane_potential_wo_spike: longblob # (mV) membrane potential without spike data, derived from membrane potential recording    
    membrane_potential_start_time: float # (s) first timepoint of membrane potential recording
    membrane_potential_sampling_rate: float # (Hz) sampling rate of membrane potential recording
    """

    def make(self, key):
        # ============ Dataset ============
        sess_data_dir = os.path.join('.', 'data', 'WholeCellData', 'Data')
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})', file=sys.stderr)
            return

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
        sess_data_dir = os.path.join('.', 'data', 'WholeCellData', 'Data')
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})', file=sys.stderr)
            return

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
class Spike(dj.Imported):
    definition = """
    -> Cell
    ---
    spike_times: longblob  # (s) time of each spike, with respect to the start of session 
    """

    def make(self, key):
        # ============ Dataset ============
        sess_data_dir = os.path.join('.', 'data', 'WholeCellData', 'Data')
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(sess_data_dir, key)
        if sess_data_file is None:
            print(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})', file=sys.stderr)
            return

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

    def make(self, key):
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get raw
        fs, first_time_point, Vm_wo_spike, Vm_w_spike = (MembranePotential & key).fetch1(
            'membrane_potential_sampling_rate', 'membrane_potential_start_time', 'membrane_potential_wo_spike',
            'membrane_potential')
        # segmentation
        segmented_Vm_wo_spike = analysis.perform_trial_segmentation(key, event_name, pre_stim_dur, post_stim_dur,
                                                           Vm_wo_spike, fs, first_time_point)
        segmented_Vm_w_spike = analysis.perform_trial_segmentation(key, event_name, pre_stim_dur, post_stim_dur,
                                                          Vm_w_spike, fs, first_time_point)
        self.insert1(dict(key,
                          segmented_mp = segmented_Vm_w_spike,
                          segmented_mp_wo_spike = segmented_Vm_wo_spike))
        print(f'Perform trial-segmentation of membrane potential for trial: {key["trial_id"]}')


@schema
class TrialSegmentedCurrentInjection(dj.Computed):
    definition = """
    -> CurrentInjection
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_current_injection: longblob
    """

    def make(self, key):
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        # get raw
        fs, first_time_point, current_injection = (CurrentInjection & key).fetch1(
            'current_injection_sampling_rate', 'current_injection_start_time', 'current_injection')
        segmented_current_injection = analysis.perform_trial_segmentation(key, event_name, pre_stim_dur, post_stim_dur,
                                                                 current_injection, fs, first_time_point)
        self.insert1(dict(key, segmented_current_injection = segmented_current_injection))
        print(f'Perform trial-segmentation of current injection for trial: {key["trial_id"]}')
