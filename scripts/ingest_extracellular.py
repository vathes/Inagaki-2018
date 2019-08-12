# -*- coding: utf-8 -*-

import os
import re
from datetime import datetime

import numpy as np
from decimal import Decimal
import scipy.io as sio
import pandas as pd
from tqdm import tqdm
import glob
from decimal import Decimal
import datajoint as dj

from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)
from pipeline import extracellular_path as path

# ================== Dataset ==================
# Fixex-delay
fixed_delay_xlsx = pd.read_excel(
    os.path.join(path, 'FixedDelayTask', 'SI_table_2_bilateral_perturb.xlsx'),
    index_col =0, usecols='A, P, Q, R, S', skiprows=2, nrows=20)
fixed_delay_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'session_time']
fixed_delay_xlsx['sex'] = 'Unknown'
fixed_delay_xlsx['sess_type'] = 'Auditory task'
fixed_delay_xlsx['delay_duration'] = 2
# Random-long-delay
random_long_delay_xlsx = pd.read_excel(
    os.path.join(path, 'RandomDelayTask', 'SI_table_3_random_delay_perturb.xlsx'),
    index_col =0, usecols='A, P, Q, R, S', skiprows=5, nrows=23)
random_long_delay_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'session_time']
random_long_delay_xlsx['sex'] = 'Unknown'
random_long_delay_xlsx['sess_type'] = 'Auditory task'
random_long_delay_xlsx['delay_duration'] = np.nan
# Random-short-delay
random_short_delay_xlsx = pd.read_excel(
    os.path.join(path, 'RandomDelayTask', 'SI_table_3_random_delay_perturb.xlsx'),
    index_col =0, usecols='A, F, G, H, I', skiprows=42, nrows=11)
random_short_delay_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'session_time']
random_short_delay_xlsx['sex'] = 'Unknown'
random_short_delay_xlsx['sess_type'] = 'Auditory task'
random_short_delay_xlsx['delay_duration'] = np.nan
# Tactile-task
tactile_xlsx = pd.read_csv(
    os.path.join(path, 'TactileTask', 'Whisker_taskTavle_for_paper.csv'),
    index_col =0, usecols= [0, 5, 6, 7, 8, 9], skiprows=1, nrows=30)
tactile_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'sex', 'session_time']
tactile_xlsx = tactile_xlsx.reindex(columns=['subject_id', 'genotype', 'date_of_birth', 'session_time', 'sex'])
tactile_xlsx['sess_type'] = 'Tactile task'
tactile_xlsx['delay_duration'] = 1.2
# Sound-task 1.2s
sound12_xlsx = pd.read_csv(
    os.path.join(path, 'Sound task 1.2s', 'OppositeTask12_for_paper.csv'),
    index_col =0, usecols= [0, 5, 6, 7, 8, 9], skiprows=1, nrows=37)
sound12_xlsx.columns = ['subject_id', 'genotype', 'date_of_birth', 'sex', 'session_time']
sound12_xlsx = sound12_xlsx.reindex(columns=['subject_id', 'genotype', 'date_of_birth', 'session_time', 'sex'])
sound12_xlsx['sess_type'] = 'Auditory task'
sound12_xlsx['delay_duration'] = 1.2
# concat all 5
meta_data = pd.concat([fixed_delay_xlsx, random_long_delay_xlsx, random_short_delay_xlsx, tactile_xlsx, sound12_xlsx])

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
fnames = np.hstack(glob.glob(os.path.join(dir_files[0], '*.mat'))
                   for dir_files in os.walk(path) if len(dir_files[1]) == 0)

for fname in fnames:
    mat = sio.loadmat(fname, struct_as_record=False, squeeze_me=True)
    mat_units = mat['unit']
    mat_trial_info = mat.get('trial_info')
    this_sess = meta_data.loc[re.sub('_units.mat|_JRC_units', '', os.path.split(fname)[-1])]
    print(f'\nReading: {this_sess.name}')

    subject_info = dict(subject_id=this_sess.subject_id.lower(),
                        date_of_birth=utilities.parse_date(str(this_sess.date_of_birth)),
                        sex=this_sess.sex[0].upper(),
                        species='Mus musculus',  # not available, hard-coded here
                        animal_source='N/A')  # animal source not available from data

    allele_dict = {alias.lower(): allele for alias, allele in subject.AlleleAlias.fetch()}
    regex_str = '|'.join([re.escape(alias) for alias in allele_dict.keys()])
    alleles = [allele_dict[s.lower()] for s in re.findall(regex_str, this_sess.genotype, re.I)]

    with subject.Subject.connection.transaction:
        if subject_info not in subject.Subject.proj():
            subject.Subject.insert1(subject_info, ignore_extra_fields=True)
            subject.Subject.Allele.insert((dict(subject_info, allele=k)
                                           for k in alleles), ignore_extra_fields = True)

    # ==================== session ====================
    # -- session_time
    session_time = utilities.parse_date(str(this_sess.session_time))
    session_info = dict(subject_info,
                        session_id='_'.join(this_sess.name.split('_')[:2]),
                        session_time=session_time)

    experimenters = ['Hidehiko Inagaki']  # hard-coded here
    experiment_types = this_sess.sess_type
    experiment_types = [experiment_types] if isinstance(experiment_types, str) else experiment_types
    experiment_types.append('extracellular')

    # experimenter and experiment type (possible multiple experimenters or types)
    # no experimenter info
    acquisition.ExperimentType.insert(zip(experiment_types), skip_duplicates=True)

    with acquisition.Session.connection.transaction:
        if session_info not in acquisition.Session.proj():
            acquisition.Session.insert1(session_info, ignore_extra_fields=True)
            acquisition.Session.Experimenter.insert((dict(session_info, experimenter=k)
                                                     for k in experimenters),
                                                    ignore_extra_fields=True)
            acquisition.Session.ExperimentType.insert((dict(session_info, experiment_type=k)
                                                       for k in experiment_types),
                                                      ignore_extra_fields=True)
        print(f'\nCreating Session - Subject: {subject_info["subject_id"]} - Date: {session_info["session_time"]}')

    # ==================== Trials ====================
    # Trial Info for all units are the same -> pick unit[0] to extract trial info
    unit_0 = mat_units[0]
    trial_key = dict(session_info, trial_counts=len(unit_0.Trial_info.Trial_types))

    with acquisition.TrialSet.connection.transaction:
        if trial_key not in acquisition.TrialSet.proj():
            fs = unit_0.Meta_data.parameters.Sample_Rate
            # handle different fieldnames "Sampling_start" vs "Sample_start"
            if 'Sample_start' not in unit_0.Behavior._fieldnames and 'Sampling_start' in unit_0.Behavior._fieldnames:
                unit_0.Behavior.Sample_start = unit_0.Behavior.Sampling_start
            if unit_0.Behavior.stim_trial_vector.size == 0:
                unit_0.Behavior.stim_trial_vector = [True if re.search('_s_', str(tr_type)) else False
                                                     for tr_type in unit_0.Trial_info.Trial_types]
            # compute delay_duration
            delay_dur = np.nanmedian(unit_0.Behavior.Cue_start - unit_0.Behavior.Delay_start)

            print('\nInsert trial information')
            acquisition.TrialSet.insert1(trial_key, allow_direct_insert=True, ignore_extra_fields = True)

            for tr_idx, (stim_trial, trial_type_of_response, trial_type,
                         first_lick, cue_start, delay_start, sample_start) in tqdm(
                enumerate(zip(unit_0.Behavior.stim_trial_vector, unit_0.Behavior.Trial_types_of_response_vector,
                              unit_0.Trial_info.Trial_types, unit_0.Behavior.First_lick, unit_0.Behavior.Cue_start,
                              unit_0.Behavior.Delay_start, unit_0.Behavior.Sample_start))):

                trial_key['trial_id'] = tr_idx + 1  # trial-number starts from 1
                trial_key['start_time'] = mat_trial_info[tr_idx].onset / fs if mat_trial_info is not None else None  # hard-coded here, no trial-start times found in data for 1st paper
                trial_key['stop_time'] = mat_trial_info[tr_idx].offset / fs if mat_trial_info is not None else None  # hard-coded here, no trial-end times found in data
                trial_key['trial_stim_present'] = bool(stim_trial != 0)
                trial_key['trial_is_good'] = bool(unit_0.Trial_info.Trial_range_to_analyze[0]
                                                  <= tr_idx <= unit_0.Trial_info.Trial_range_to_analyze[-1])
                trial_key['trial_type'], trial_key['trial_response'] = trial_type_and_response_dict[trial_type_of_response]
                trial_key['delay_duration'] = Decimal(cue_start - delay_start).quantize(Decimal('0.1'))
                acquisition.TrialSet.Trial.insert1(trial_key, ignore_extra_fields = True, skip_duplicates = True,
                                                   allow_direct_insert = True)

                # ======== Now add trial event timing to the EventTime part table ====
                events_time = dict(trial_start=0,
                                   trial_stop=(trial_key['stop_time'] - trial_key['start_time']
                                               if mat_trial_info is not None else None),
                                   first_lick=first_lick,
                                   cue_start=cue_start,
                                   delay_start=delay_start,
                                   sampling_start=sample_start)
                # -- events timing
                acquisition.TrialSet.EventTime.insert((dict(trial_key, trial_event=k, event_time=e)
                                                       for k, e in events_time.items()),
                                                      ignore_extra_fields = True, skip_duplicates = True,
                                                      allow_direct_insert = True)
                # ======== Now add trial stimulation descriptors to the TrialPhotoStimInfo table ====
                trial_key['photo_stim_period'] = 'early delay'  # TODO: hardcoded here because this info is not available from data
                trial_key['photo_stim_power'] = (re.search('(?<=_)\d+(?=mW_)', str(trial_type)).group()  # str() to safeguard against np.array([]) (probably typo)
                                                 if re.search('(?<=_)\d+(?=mW_)', str(trial_type)) else None)
                stimulation.TrialPhotoStimParam.insert1(trial_key, ignore_extra_fields=True, allow_direct_insert=True)

    # ==================== Extracellular ====================
    # no info about Probe or recording location from data, all hardcoded from paper
    channel_counts = 64
    chn_per_shank = 32
    probe_name = 'A2x32-8mm-25-250-165'
    # -- Probe
    with reference.Probe.connection.transaction:
        if {'probe_name': probe_name, 'channel_counts': channel_counts} not in reference.Probe.proj():
            reference.Probe.insert1({'probe_name': probe_name, 'channel_counts': channel_counts})
            reference.Probe.Channel.insert({'probe_name': probe_name, 'channel_counts': channel_counts,
                                            'channel_id': ch_idx, 'shank_id': int(ch_idx <= chn_per_shank) + 1}
                                           for ch_idx in np.arange(channel_counts) + 1)

    brain_region = 'ALM'
    hemisphere = 'left'
    brain_location = {'brain_region': brain_region,
                      'brain_subregion': 'N/A',
                      'cortical_layer': 'N/A',
                      'hemisphere': hemisphere}
    reference.BrainLocation.insert1(brain_location, skip_duplicates=True)

    # -- ProbeInsertion
    probe_insertion = dict({**session_info, **brain_location},
                           probe_name=probe_name, channel_counts=channel_counts)
    extracellular.ProbeInsertion.insert1(probe_insertion, ignore_extra_fields=True, skip_duplicates=True)

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
                                device_name = stim_device,
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
    if dict(session_info, photostim_datetime=session_info['session_time']) not in stimulation.PhotoStimulation.proj():
        stimulation.PhotoStimulation.insert1(dict({**session_info, **action_location, **photim_stim_protocol},
                                                  photostim_datetime=session_info['session_time'])
                                             , ignore_extra_fields=True)


