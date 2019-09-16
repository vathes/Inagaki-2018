# -*- coding: utf-8 -*-

import os
import re
from datetime import datetime

import numpy as np
from decimal import Decimal
import scipy.io as sio
import pandas as pd
from tqdm import tqdm
import uuid
import datajoint as dj
import glob

from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)
from pipeline import intracellular_path as path

# ================== Dataset ==================
xlsname = 'SI_table_1_wc_cell_list.xlsx'
meta_data = pd.read_excel(os.path.join(path, xlsname),
                          index_col =0,
                          usecols='A:N',
                          skiprows=[1],
                          nrows=79)
meta_data.columns = ['experiment_type', 'depth_um',
                     'correct_contra_trial_count', 'correct_ipsi_trial_count',
                     'performance',
                     'delay_selectivity_SR', 'delay_selectivity_Vm',
                     'subject_id', 'genotype', 'date_of_birth', 'session_time',
                     'current_injection', 'from_guo_inagaki_2017']

trial_type_and_response_dict = {1: ('lick right', 'correct'),
                                2: ('lick left', 'correct'),
                                3: ('lick right', 'incorrect'),
                                4: ('lick left', 'incorrect'),
                                5: ('lick right', 'early lick'),
                                6: ('lick left', 'early lick'),
                                7: ('lick right', 'early lick'),
                                8: ('lick left', 'early lick'),
                                9: ('lick right', 'no response'),
                                10: ('lick left', 'no response'),
                                11: ('N/A', 'N/A')}

# ========================== METADATA ==========================
# ==================== subject ====================
for fname in glob.glob(os.path.join(path, 'Data', '*.mat')):
    mat_data = sio.loadmat(fname, struct_as_record = False, squeeze_me = True)['wholeCell']
    fname = (os.path.split(fname)[-1]).replace('.mat', '')
    this_sess = meta_data.loc[f'Cell {mat_data.cell_id}']
    print(f'\nReading: {fname}')

    subject_info = dict(subject_id=this_sess.subject_id.lower(),
                        date_of_birth=utilities.parse_date(str(this_sess.date_of_birth)),
                        species='Mus musculus',  # not available, hard-coded here
                        animal_source='N/A')  # animal source not available from data, nor 'sex'

    # allele
    allele_dict = {alias.lower(): allele for alias, allele in subject.AlleleAlias.fetch()}
    regex_str = '|'.join([re.escape(alias) for alias in allele_dict.keys()])
    alleles = [allele_dict[s.lower()] for s in re.findall(regex_str, this_sess.genotype, re.I)]

    with subject.Subject.connection.transaction:
        if subject_info not in subject.Subject.proj():
            subject.Subject.insert1(subject_info, ignore_extra_fields=True)
            subject.Subject.Allele.insert((dict(subject_info, allele = k)
                                           for k in alleles), ignore_extra_fields = True)

    # ==================== session ====================
    # -- session_time
    session_time = utilities.parse_date(str(this_sess.session_time))
    session_info = dict(subject_info,
                        session_id=fname,
                        session_time=session_time)

    experimenters = ['Hidehiko Inagaki']  # hard-coded here
    experiment_types = this_sess.experiment_type
    experiment_types = [experiment_types] if isinstance(experiment_types, str) else experiment_types
    experiment_types.append('intracellular')
    experiment_types.append(re.search("regular|EPSP", fname).group())


    # experimenter and experiment type (possible multiple experimenters or types)
    # no experimenter info
    acquisition.ExperimentType.insert(zip(experiment_types), skip_duplicates=True)

    with acquisition.Session.connection.transaction:
        if session_info not in acquisition.Session.proj():
            acquisition.Session.insert1(session_info, ignore_extra_fields=True)
            acquisition.Session.Experimenter.insert((dict(session_info, experimenter=k)
                                                     for k in experimenters), ignore_extra_fields=True)
            acquisition.Session.ExperimentType.insert((dict(session_info, experiment_type=k)
                                                       for k in experiment_types), ignore_extra_fields=True)
        print(f'\nCreating Session - Subject: {subject_info["subject_id"]} - Date: {session_info["session_time"]}')

    # ==================== Intracellular ====================
    # no info on recording device, brain location of cell available from data, hard-coded from the paper
    ie_device = 'Multiclamp 700B'
    brain_region = 'ALM'
    hemisphere = 'left'
    brain_location = {'brain_region': brain_region,
                      'brain_subregion': 'N/A',
                      'cortical_layer': 'N/A',
                      'hemisphere': hemisphere}
    reference.BrainLocation.insert1(brain_location, skip_duplicates=True)

    # -- Cell
    cell_id = f'cell_{mat_data.cell_id}'
    cell_key = dict({**session_info, **brain_location},
                    cell_id=cell_id,
                    cell_type=mat_data.Pyr_or_GABA,
                    cell_depth=float(this_sess.depth_um),
                    device_name=ie_device)
    intracellular.Cell.insert1(cell_key, ignore_extra_fields=True, skip_duplicates=True)

    # ==================== Trials ====================
    fs = mat_data.recording_data.sample_rate
    trial_key = dict(session_info, trial_counts=len(mat_data.behavioral_data.behav_timing))

    with acquisition.TrialSet.connection.transaction:
        if trial_key not in acquisition.TrialSet.proj():
            print('\nInsert trial information')
            acquisition.TrialSet.insert1(trial_key, allow_direct_insert = True, ignore_extra_fields = True)

            for tr_idx, events_time in tqdm(enumerate(mat_data.behavioral_data.behav_timing)):
                trial_key['trial_id'] = tr_idx + 1  # trial-number starts from 1
                trial_key['start_time'] = mat_data.behavioral_data.trial_onset_bin[tr_idx] / fs
                trial_key['stop_time'] = min(trial_key['start_time'] + events_time.end_time,
                                             (len(mat_data.recording_data.Vm) - 1) / fs)
                trial_key['trial_stim_present'] = bool(mat_data.behavioral_data.AOM_on_or_off[tr_idx])
                trial_key['trial_is_good'] = True  # no info of trial good/bad status, assuming all trials are good
                trial_key['trial_type'], trial_key['trial_response'] = trial_type_and_response_dict[
                    mat_data.behavioral_data.trial_type_vector[tr_idx]]
                trial_key['delay_duration'] = 1.2  # hard-coded here (the same for whole cell)
                acquisition.TrialSet.Trial.insert1(trial_key, ignore_extra_fields = True, skip_duplicates = True,
                                                   allow_direct_insert = True)

                # ======== Now add trial event timing to the EventTime part table ====
                events = dict(
                    events_time.__dict__,
                    trial_start = 0,
                    trial_stop = trial_key['stop_time'] - trial_key['start_time'],
                    first_lick = min([np.array(getattr(events_time, l)).flatten()[
                                          0]  # events_time(l) could be empty ([]), a single time (float) or multiple times (array)
                                      if np.array(getattr(events_time, l)).flatten().size > 0 else np.nan
                                      for l in ('lickL_on_time', 'lickR_on_time')]),
                    current_injection_start = mat_data.behavioral_data.tail_current_injection_onset_bin[tr_idx] / fs)
                # -- events timing
                acquisition.TrialSet.EventTime.insert((dict(trial_key, trial_event = k, event_time = events[k])
                                                       for k in ['trial_start', 'trial_stop', 'cue_start',
                                                                 'cue_end', 'sampling_start', 'delay_start',
                                                                 'current_injection_start', 'first_lick']),
                                                      ignore_extra_fields = True, skip_duplicates = True,
                                                      allow_direct_insert = True)

    # ==================== photostim ====================
    # no info on photostim available from data, all photostim info here are hard-coded from the paper
    brain_region = 'ALM'
    hemisphere = 'bilateral'
    coord_ap_ml_dv = [2.5, 1.5, 0]
    stim_device = 'laser'  # hard-coded here..., could not find a more specific name from metadata
    device_desc = 'laser (Laser Quantum) were controlled by an acousto-optical modulator (AOM; Quanta Tech) and a shutter (Vincent Associates)'

    # -- Device
    stimulation.PhotoStimDevice.insert1({'device_name': stim_device, 'device_desc': device_desc}, skip_duplicates=True)

    photim_stim_protocol = dict(protocol=1,
                                device_name=stim_device,
                                photo_stim_excitation_lambda=473,
                                photo_stim_notes='photostimulate four spots in each hemisphere, centered on ALM (AP 2.5 mm; ML 1.5 mm) with 1 mm spacing (in total eight spots bilaterally) using scanning Galvo mirrors',
                                photo_stim_duration=1000,
                                photo_stim_freq=40,
                                photo_stim_shape='sinusoidal',
                                photo_stim_method='laser')

    # -- BrainLocation
    brain_location = {'brain_region': brain_region,
                      'brain_subregion': 'N/A',
                      'cortical_layer': 'N/A',
                      'hemisphere': hemisphere}
    reference.BrainLocation.insert1(brain_location, skip_duplicates=True)
    # -- ActionLocation
    action_location = dict(brain_location,
                           coordinate_ref = 'bregma',
                           coordinate_ap = round(Decimal(coord_ap_ml_dv[0]), 2),
                           coordinate_ml = round(Decimal(coord_ap_ml_dv[1]), 2),
                           coordinate_dv = round(Decimal(coord_ap_ml_dv[2]), 2))
    reference.ActionLocation.insert1(action_location, skip_duplicates=True)

    # -- PhotoStimulationProtocol
    stimulation.PhotoStimProtocol.insert1(photim_stim_protocol, skip_duplicates=True)
    # -- PhotoStimulation
    # only 1 photostim per session, perform at the same time with session
    photo_stimulation = dict({**session_info, **photim_stim_protocol}, photostim_datetime=session_info['session_time'])
    if photo_stimulation not in stimulation.PhotoStimulation.proj():
        stimulation.PhotoStimulation.insert1({**photo_stimulation,
                                              **action_location,
                                              **(dict(photostim_start_time=0,
                                                      photostim_timeseries=mat_data.recording_data.AOM,
                                                      photostim_sampling_rate=mat_data.recording_data.sample_rate)
                                                 if mat_data.recording_data.AOM.size > 0 else dict())},
                                             ignore_extra_fields=True)
