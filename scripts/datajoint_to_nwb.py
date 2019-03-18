#!/usr/bin/env python3
import os

import sys
from datetime import datetime
from dateutil.tz import tzlocal
import pytz
import re
import numpy as np
import pandas as pd
import warnings
import tqdm

from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)
import pynwb
from pynwb import NWBFile, NWBHDF5IO

warnings.filterwarnings('ignore', module='pynwb')
# =============================================
# Each NWBFile represent a session, thus for every session in acquisition.Session, we build one NWBFile

for session_key in tqdm.tqdm(acquisition.Session.fetch('KEY')):
    this_session = (acquisition.Session & session_key).fetch1()
    # =============== General ====================
    # -- NWB file - a NWB2.0 file for each session
    nwbfile = NWBFile(
        session_description=this_session['session_note'],
        identifier='_'.join([this_session['subject_id'],
                             this_session['session_time'].strftime('%Y-%m-%d'),
                             this_session['session_id']]),
        session_start_time=this_session['session_time'],
        file_create_date=datetime.now(tzlocal()),
        experimenter='; '.join((acquisition.Session.Experimenter & session_key).fetch('experimenter')),
        institution='Janelia Research Campus',  # TODO: not in pipeline
        related_publications='')  # TODO: not in pipeline
    # -- subject
    subj = (subject.Subject & session_key).fetch1()
    nwbfile.subject = pynwb.file.Subject(
        subject_id=this_session['subject_id'],
        description=subj['subject_description'],
        genotype=' x '.join((subject.Subject.Allele & session_key).fetch('allele')),
        sex=subj['sex'],
        species=subj['species'])
    # =============== Intracellular ====================
    cell = ((intracellular.Cell & session_key).fetch1()
            if len(intracellular.Cell & session_key) == 1
            else None)
    if cell:
        # metadata
        whole_cell_device = nwbfile.create_device(name=cell['device_name'])
        ic_electrode = nwbfile.create_ic_electrode(
            name=cell['cell_id'],
            device=whole_cell_device,
            description='N/A',
            filtering='N/A',  # TODO: not in pipeline
            location='; '.join([f'{k}: {str(v)}'
                                for k, v in dict((reference.BrainLocation & cell).fetch1(),
                                                 depth=cell['cell_depth']).items()]))
        # acquisition - membrane potential
        mp, mp_wo_spike, mp_start_time, mp_fs = (intracellular.MembranePotential & cell).fetch1(
            'membrane_potential', 'membrane_potential_wo_spike',
            'membrane_potential_start_time', 'membrane_potential_sampling_rate')
        nwbfile.add_acquisition(pynwb.icephys.PatchClampSeries(name='membrane_potential',
                                                               electrode=ic_electrode,
                                                               unit='mV',  # TODO: not in pipeline
                                                               conversion=1e-3,
                                                               gain=1.0,  # TODO: not in pipeline
                                                               data=mp,
                                                               starting_time=mp_start_time,
                                                               rate=mp_fs))
        # acquisition - current injection
        if (intracellular.CurrentInjection & cell):
            current_injection, ci_start_time, ci_fs = (intracellular.CurrentInjection & cell).fetch1(
                'current_injection', 'current_injection_start_time', 'current_injection_sampling_rate')
            nwbfile.add_stimulus(pynwb.icephys.CurrentClampStimulusSeries(name='current_injection',
                                                                          electrode=ic_electrode,
                                                                          unit='nA',
                                                                          conversion=1e-6,
                                                                          gain=1.0,
                                                                          data=current_injection,
                                                                          starting_time=ci_start_time,
                                                                          rate=ci_fs))

        # analysis - membrane potential without spike
        mp_rmv_spike = nwbfile.create_processing_module(name='membrane_potential_spike_removal',
                                                        description='Spike removal')
        mp_rmv_spike.add_data_interface(pynwb.icephys.PatchClampSeries(name='membrane_potential_without_spike',
                                                                       electrode=ic_electrode,
                                                                       unit='mV',
                                                                       conversion=1e-3,
                                                                       gain=1.0,
                                                                       data=mp_wo_spike,
                                                                       starting_time=mp_start_time,
                                                                       rate=mp_fs))

    # =============== Extracellular ====================
    probe_insertion = ((extracellular.ProbeInsertion & session_key).fetch1()
                       if extracellular.ProbeInsertion & session_key
                       else None)
    if probe_insertion:
        probe = nwbfile.create_device(name = probe_insertion['probe_name'])
        electrode_group = nwbfile.create_electrode_group(
            name='; '.join([f'{probe_insertion["probe_name"]}: {str(probe_insertion["channel_counts"])}']),
            description = 'N/A',
            device = probe,
            location = '; '.join([f'{k}: {str(v)}' for k, v in
                                  (reference.BrainLocation & probe_insertion).fetch1().items()]))

        for chn in (reference.Probe.Channel & probe_insertion).fetch(as_dict=True):
            nwbfile.add_electrode(id=chn['channel_id'],
                                  group=electrode_group,
                                  filtering='',  # TODO: not in pipeline
                                  imp=-1.,  # TODO: not in pipeline
                                  x=0.0,  # not available from data
                                  y=0.0,  # not available from data
                                  z=0.0,  # not available from data
                                  location=electrode_group.location)

        # --- unit spike times ---
        nwbfile.add_unit_column(name='depth', description='depth this unit')
        nwbfile.add_unit_column(name='spike_width', description='spike width of this unit')
        nwbfile.add_unit_column(name='cell_type', description='cell type (e.g. wide width, narrow width spiking)')

        for unit in (extracellular.UnitSpikeTimes & probe_insertion).fetch(as_dict=True):
            # make an electrode table region (which electrode(s) is this unit coming from)
            nwbfile.add_unit(id=unit['unit_id'],
                             electrodes=(unit['channel_id']
                                         if isinstance(unit['channel_id'], np.ndarray) else [unit['channel_id']]),
                             depth=unit['unit_depth'],
                             spike_width=unit['unit_spike_width'],
                             cell_type=unit['unit_cell_type'],
                             spike_times=unit['spike_times'],
                             waveform_mean=unit['spike_waveform'])

    # =============== Behavior ====================
    # Note: for this study, raw behavioral data were not available, only trialized data were provided
    # here, we reconstruct raw behavioral data by concatenation
    trial_seg_setting = (analysis.TrialSegmentationSetting & 'trial_seg_setting=0').fetch1()
    seg_behav_query = (behavior.TrialSegmentedLickTrace * acquisition.TrialSet.Trial
                 * (analysis.RealignedEvent.RealignedEventTime & 'trial_event="trial_start"')
                 & session_key & trial_seg_setting)

    if seg_behav_query:
        behav_acq = pynwb.behavior.BehavioralTimeSeries(name = 'lick_times')
        nwbfile.add_acquisition(behav_acq)
        seg_behav = pd.DataFrame(seg_behav_query.fetch('start_time', 'realigned_event_time',
                                                       'segmented_lick_left_on',
                                                       'segmented_lick_left_off',
                                                       'segmented_lick_right_on',
                                                       'segmented_lick_right_off')).T
        seg_behav.columns = ['start_time', 'realigned_event_time', 'segmented_lick_left_on',
                             'segmented_lick_left_off', 'segmented_lick_right_on', 'segmented_lick_right_off']
        for behav_name in ['lick_left_on', 'lick_left_off', 'lick_right_on', 'lick_right_off']:
            lick_times = np.hstack(r['segmented_'+behav_name] - r.realigned_event_time + r.start_time
                                 for _, r in seg_behav.iterrows())
            behav_acq.create_timeseries(
                name = behav_name,
                unit = 'a.u.',
                conversion = 1.0,
                data = np.full_like(lick_times, 1),
                timestamps = lick_times)

    # =============== Photostimulation ====================
    photostim = ((stimulation.PhotoStimulation & session_key).fetch1()
                       if stimulation.PhotoStimulation & session_key
                       else None)
    if photostim:
        photostim_device = (stimulation.PhotoStimDevice & photostim).fetch1()
        stim_device = nwbfile.create_device(name=photostim_device['device_name'])
        stim_site = pynwb.ogen.OptogeneticStimulusSite(
            name='-'.join([photostim['hemisphere'], photostim['brain_region']]),
            device=stim_device,
            excitation_lambda=float((stimulation.PhotoStimProtocol & photostim).fetch1('photo_stim_excitation_lambda')),
            location = '; '.join([f'{k}: {str(v)}' for k, v in
                                  (reference.ActionLocation & photostim).fetch1().items()]),
            description=(stimulation.PhotoStimProtocol & photostim).fetch1('photo_stim_notes'))
        nwbfile.add_ogen_site(stim_site)

        if photostim['photostim_timeseries'] is not None:
            nwbfile.add_stimulus(pynwb.ogen.OptogeneticSeries(
                name='_'.join(['photostim_on', photostim['photostim_datetime'].strftime('%Y-%m-%d_%H-%M-%S')]),
                site=stim_site,
                unit = 'mW',
                resolution = 0.0,
                conversion = 1e-6,
                data = photostim['photostim_timeseries'],
                starting_time = photostim['photostim_start_time'],
                rate = photostim['photostim_sampling_rate']))

    # =============== TrialSet ====================
    # NWB 'trial' (of type dynamic table) by default comes with three mandatory attributes:
    #                                                                       'id', 'start_time' and 'stop_time'.
    # Other trial-related information needs to be added in to the trial-table as additional columns (with column name
    # and column description)
    if acquisition.TrialSet & session_key:
        # Get trial descriptors from TrialSet.Trial and TrialStimInfo
        trial_columns = [{'name': tag,
                          'description': re.sub('\s+:|\s+', ' ', re.search(
                              f'(?<={tag})(.*)', str((acquisition.TrialSet.Trial * stimulation.TrialPhotoStimParam).heading)).group())}
                         for tag in (acquisition.TrialSet.Trial * stimulation.TrialPhotoStimParam).fetch(as_dict=True, limit=1)[0].keys()
                         if tag not in (acquisition.TrialSet.Trial & stimulation.TrialPhotoStimParam).primary_key + ['start_time', 'stop_time']]

        # Trial Events - discard 'trial_start' and 'trial_stop' as we already have start_time and stop_time
        trial_events = set(((acquisition.TrialSet.EventTime & session_key)
                            - [{'trial_event': 'trial_start'}, {'trial_event': 'trial_stop'}]).fetch('trial_event'))
        event_names = [{'name': e, 'description': d}
                       for e, d in zip(*(reference.ExperimentalEvent & [{'event': k}
                                                                        for k in trial_events]).fetch('event',
                                                                                                      'description'))]
        # Add new table columns to nwb trial-table for trial-label
        for c in trial_columns + event_names:
            nwbfile.add_trial_column(**c)

        photostim_tag_default = {tag: '' for tag in stimulation.TrialPhotoStimParam().fetch(as_dict=True, limit=1)[0].keys()
                                 if tag not in stimulation.TrialPhotoStimParam.primary_key}
        # Add entry to the trial-table
        for trial in (acquisition.TrialSet.Trial & session_key).fetch(as_dict=True):
            events = dict(zip(*(acquisition.TrialSet.EventTime & trial
                                & [{'trial_event': e} for e in trial_events]).fetch('trial_event', 'event_time')))
            photostim_tag = (stimulation.TrialPhotoStimParam & trial).fetch(as_dict=True)
            trial_tag_value = ({**trial, **events, **photostim_tag[0]}
                               if len(photostim_tag) == 1 else {**trial, **events, **photostim_tag_default})

            trial_tag_value['id'] = trial_tag_value['trial_id']  # rename 'trial_id' to 'id'
            trial_tag_value['delay_duration'] = float(trial_tag_value['delay_duration'])  # convert Decimal to float
            # convert None to np.nan since nwb fields does not take None
            for k, v in trial_tag_value.items():
                trial_tag_value[k] = v if v is not None else np.nan
            [trial_tag_value.pop(k) for k in acquisition.TrialSet.Trial.primary_key]
            nwbfile.add_trial(**trial_tag_value)

    # =============== Write NWB 2.0 file ===============
    save_path = os.path.join('data', 'NWB 2.0')
    save_file_name = ''.join([nwbfile.identifier, '.nwb'])
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    with NWBHDF5IO(os.path.join(save_path, save_file_name), mode = 'w') as io:
        io.write(nwbfile)
        print(f'Write NWB 2.0 file: {save_file_name}')
        
