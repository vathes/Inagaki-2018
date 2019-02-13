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

from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)

# ================== Dataset ==================
path = os.path.join('.', 'data', 'SiliconProbeData')

# Fixex-delay
fixed_delay_xlsx = pd.read_excel(
    os.path.join('.', 'data', 'SiliconProbeData', 'FixedDelayTask', 'SI_table_2_bilateral_perturb.xlsx'),
    index_col =0, usecols='A, P, Q, R, S', skiprows=2, nrows=20)
fixed_delay_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'session_time']
fixed_delay_xlsx['sess_type'] = 'fixed_delay'
# Random-long-delay
random_long_delay_xlsx = pd.read_excel(
    os.path.join('.', 'data', 'SiliconProbeData', 'RandomDelayTask', 'SI_table_3_random_delay_perturb.xlsx'),
    index_col =0, usecols='A, P, Q, R, S', skiprows=5, nrows=23)
random_long_delay_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'session_time']
random_long_delay_xlsx['sess_type'] = 'random_long_delay'
# Random-short-delay
random_short_delay_xlsx = pd.read_excel(
    os.path.join('.', 'data', 'SiliconProbeData', 'RandomDelayTask', 'SI_table_3_random_delay_perturb.xlsx'),
    index_col =0, usecols='A, F, G, H, I', skiprows=42, nrows=11)
random_short_delay_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'session_time']
random_short_delay_xlsx['sess_type'] = 'random_short_delay'
# concat all 3
meta_data = pd.concat([fixed_delay_xlsx, random_long_delay_xlsx, random_short_delay_xlsx])

trial_type_and_response_dict = {1: ('lick right', 'correct'),
                                2: ('lick left', 'correct'),
                                3: ('lick right', 'incorrect'),
                                4: ('lick left', 'incorrect'),
                                5: ('lick right', 'no response'),
                                6: ('lick left', 'no response'),
                                7: ('lick right', 'early lick'),
                                8: ('lick left', 'early lick'),
                                9: ('lick right', 'early lick'),
                                10: ('lick left', 'early lick'),
                                0: ('photo-tagging', 'N/A')}

# ========================== METADATA ==========================
# ==================== subject ====================
fnames = np.hstack([os.path.join(dir_files[0], f) for f in dir_files[2] if f.find('.mat') != -1]
          for dir_files in os.walk(path) if len(dir_files[1]) == 0)

for fname in fnames:
    mat_units = sio.loadmat(fname, struct_as_record = False, squeeze_me = True)['unit']
    this_sess = meta_data.loc[os.path.split(fname)[-1].replace('_units.mat', '')]
    print(f'\nReading: {this_sess.name}')

    subject_info = dict(subject_id=this_sess.subject_id.lower(),
                        date_of_birth=datetime.strptime(str(this_sess.date_of_birth), '%Y%m%d'),
                        species='Mus musculus',  # not available, hard-coded here
                        animal_source='N/A')  # animal source not available from data, nor 'sex'

    strain_dict = {alias.lower(): strain for alias, strain in subject.StrainAlias.fetch()}
    regex_str = ''.join([re.escape(alias) + '|' for alias in strain_dict.keys()])[:-1]
    strains = [strain_dict[s.lower()] for s in re.findall(regex_str, this_sess.genotype, re.I)]

    if subject_info not in subject.Subject.proj():
        with subject.Subject.connection.transaction:
            subject.Subject.insert1(subject_info, ignore_extra_fields=True)
            subject.Subject.Strain.insert((dict(subject_info, strain = k)
                                           for k in strains), ignore_extra_fields = True)

    # ==================== session ====================
    # -- session_time
    session_time = datetime.strptime(str(this_sess.session_time), '%Y%m%d')
    session_info = dict(subject_info,
                        session_id=this_sess.name,
                        session_time=session_time)

    experimenters = ['Hidehiko Inagaki']  # hard-coded here
    experiment_types = this_sess.sess_type
    experiment_types = [experiment_types] if isinstance(experiment_types, str) else experiment_types
    experiment_types.append('extracellular')

    # experimenter and experiment type (possible multiple experimenters or types)
    # no experimenter info
    acquisition.ExperimentType.insert(zip(experiment_types), skip_duplicates=True)

    if session_info not in acquisition.Session.proj():
        with acquisition.Session.connection.transaction:
            acquisition.Session.insert1(session_info, ignore_extra_fields=True)
            acquisition.Session.Experimenter.insert((dict(session_info, experimenter=k) for k in experimenters), ignore_extra_fields=True)
            acquisition.Session.ExperimentType.insert((dict(session_info, experiment_type=k) for k in experiment_types), ignore_extra_fields=True)
        print(f'\nCreating Session - Subject: {subject_info["subject_id"]} - Date: {session_info["session_time"]}')

    # ==================== Trials ====================
    # Trial Info for all units are the same -> pick unit[0] to extract trial info
    unit_0 = mat_units[0]
    fs = unit_0.Meta_data.parameters.Sample_Rate
    trial_key = dict(session_info, trial_counts=len(unit_0.Trial_info.Trial_types))
    if trial_key not in acquisition.TrialSet.proj():
        with acquisition.TrialSet.connection.transaction:
            print('\nInsert trial information')
            acquisition.TrialSet.insert1(trial_key, allow_direct_insert=True, ignore_extra_fields = True)

            for tr_idx in tqdm(np.arange(len(unit_0.Trial_info.Trial_types))):
                trial_key['trial_id'] = tr_idx
                trial_key['start_time'] = None  # hard-coded here, no trial-start times found in data
                trial_key['stop_time'] = None  # hard-coded here, no trial-end times found in data
                trial_key['trial_stim_present'] = bool(unit_0.Behavior.stim_trial_vector[tr_idx] != 0)
                trial_key['trial_is_good'] = bool(unit_0.Trial_info.Trial_range_to_analyze[0]
                                                  <= tr_idx <= unit_0.Trial_info.Trial_range_to_analyze[-1])
                trial_key['trial_type'], trial_key['trial_response'] = trial_type_and_response_dict[
                    unit_0.Behavior.Trial_types_of_response_vector[tr_idx]]
                trial_key['delay_duration'] = unit_0.Behavior.delay_dur[tr_idx] if 'delay_dur' in unit_0.Behavior._fieldnames else 1.2
                acquisition.TrialSet.Trial.insert1(trial_key, ignore_extra_fields = True, skip_duplicates = True,
                                                   allow_direct_insert = True)

                # ======== Now add trial event timing to the EventTime part table ====
                events_time = dict(trial_start=trial_key['start_time'],
                                   trial_stop=trial_key['stop_time'],
                                   first_lick=unit_0.Behavior.First_lick[tr_idx],
                                   cue_start=unit_0.Behavior.Cue_start[tr_idx],
                                   sampling_start=unit_0.Behavior.Delay_start[tr_idx],
                                   delay_start=unit_0.Behavior.Sample_start[tr_idx])
                # -- events timing
                acquisition.TrialSet.EventTime.insert((dict(trial_key, trial_event=k, event_time=e)
                                                       for k, e in events_time.items()),
                                                      ignore_extra_fields = True, skip_duplicates = True,
                                                      allow_direct_insert = True)
                # ======== Now add trial stimulation descriptors to the TrialPhotoStimInfo table ====
                trial_key['photo_stim_period'] = 'early delay'  # TODO: hardcoded here because this info is not available from data
                trial_key['photo_stim_power'] = (re.search('(?<=_)\d+(?=mW_)',
                                                           str(unit_0.Trial_info.Trial_types[tr_idx])).group()  # str() to safeguard against np.array([]) (probably typo)
                                                 if re.search('(?<=_)\d+(?=mW_)',
                                                              str(unit_0.Trial_info.Trial_types[tr_idx])) else None)
                stimulation.TrialPhotoStimInfo.insert1(trial_key, ignore_extra_fields=True, allow_direct_insert=True)

    # ==================== Extracellular ====================
    # no info about Probe or recording location from data, all hardcoded from paper
    channel_counts = 64
    chn_per_shank = 32
    probe_name = 'A2x32-8mm-25-250-165'
    # -- Probe
    if {'probe_name': probe_name, 'channel_counts': channel_counts} not in reference.Probe.proj():
        with reference.Probe.connection.transaction:
            reference.Probe.insert1(
                    {'probe_name': probe_name,
                     'channel_counts': channel_counts})
            reference.Probe.Channel.insert({'probe_name': probe_name, 'channel_counts': channel_counts,
                                            'channel_id': ch_idx, 'shank_id': int(ch_idx <= chn_per_shank) + 1}
                                           for ch_idx in np.arange(channel_counts) + 1)

    brain_region = 'ALM'
    hemisphere = 'left'
    brain_location = {'brain_region': brain_region,
                      'brain_subregion': 'N/A',
                      'cortical_layer': 'N/A',
                      'hemisphere': hemisphere}

    # -- ProbeInsertion
    probe_insertion = dict({**session_info, **brain_location},
                           probe_name=probe_name, channel_counts=channel_counts)
    if probe_insertion not in extracellular.ProbeInsertion.proj():
        extracellular.ProbeInsertion.insert1(probe_insertion, ignore_extra_fields=True)

    # ==================== photostim ====================
    # no info on photostim available from data, all photostim info here are hard-coded from the paper
    brain_region = 'ALM'
    hemisphere = 'bilateral'
    coord_ap_ml_dv = [2.5, 1.5, 0]
    stim_device = 'laser'  # hard-coded here..., could not find a more specific name from metadata
    device_desc = 'laser (Laser Quantum) were controlled by an acousto-optical modulator (AOM; Quanta Tech) and a shutter (Vincent Associates)'

    stim_info = dict(photo_stim_excitation_lambda=473,
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
    if brain_location not in reference.BrainLocation.proj():
        reference.BrainLocation.insert1(brain_location)
    # -- ActionLocation
    action_location = dict(brain_location,
                           coordinate_ref = 'bregma',
                           coordinate_ap = round(Decimal(coord_ap_ml_dv[0]), 2),
                           coordinate_ml = round(Decimal(coord_ap_ml_dv[1]), 2),
                           coordinate_dv = round(Decimal(coord_ap_ml_dv[2]), 2))
    if action_location not in reference.ActionLocation.proj():
        reference.ActionLocation.insert1(action_location)

    # -- Device
    if {'device_name': stim_device} not in stimulation.PhotoStimDevice.proj():
        stimulation.PhotoStimDevice.insert1({'device_name': stim_device, 'device_desc': device_desc})

    # -- PhotoStimulationInfo
    photim_stim_info = {**action_location, **stim_info, 'device_name': stim_device}
    if photim_stim_info not in stimulation.PhotoStimulationInfo.proj():
        stimulation.PhotoStimulationInfo.insert1(photim_stim_info)
    # -- PhotoStimulation
    # only 1 photostim per session, perform at the same time with session
    if dict(session_info, photostim_datetime=session_info['session_time']) not in stimulation.PhotoStimulation.proj():
        stimulation.PhotoStimulation.insert1(dict({**session_info, **photim_stim_info},
                                                  photostim_datetime=session_info['session_time'])
                                             , ignore_extra_fields=True)

# ====================== Starting import and compute procedure ======================
print('======== Populate() Routine =====')
extracellular.UnitSpikeTimes.populate()
extracellular.TrialSegmentedUnitSpikeTimes.populate()

