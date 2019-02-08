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
from pipeline import reference, subject, acquisition, stimulation, analysis #, behavior, ephys, action
from pipeline import utilities

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

    for s in subject.StrainAlias.fetch():
        if re.search(re.escape(s[0]), this_sess.genotype, re.I):
            subject_info['strain'] = (subject.StrainAlias & {'strain_alias': s[0]}).fetch1('strain')
            break

    if subject_info not in subject.Subject.proj():
        subject.Subject.insert1(subject_info, ignore_extra_fields = True)

    # ==================== session ====================
    # -- session_time
    session_time = datetime.strptime(str(this_sess.session_time), '%Y%m%d')
    session_info = dict(subject_info,
                        session_id=this_sess.name,
                        session_time=session_time)

    experimenters = ['Hidehiko Inagaki']  # hard-coded here
    experiment_types = this_sess.sess_type
    experiment_types = [experiment_types] if isinstance(experiment_types, str) else experiment_types

    # experimenter and experiment type (possible multiple experimenters or types)
    # no experimenter info
    acquisition.ExperimentType.insert(zip(experiment_types), skip_duplicates=True)

    if session_info not in acquisition.Session.proj():
        with acquisition.Session.connection.transaction:
            acquisition.Session.insert1(session_info, ignore_extra_fields=True)
            acquisition.Session.Experimenter.insert((dict(session_info, experiment_type=k) for k in experimenters), ignore_extra_fields=True)
            acquisition.Session.ExperimentType.insert((dict(session_info, experiment_type=k) for k in experiment_types), ignore_extra_fields=True)
        print(f'\nCreating Session - Subject: {subject_info["subject_id"]} - Date: {session_info["session_time"]}')

    # ==================== Extracellular ====================


    # ==================== Trials ====================
    # Trial Info for all units are the same -> pick unit[0] to extract trial info
    unit_0 = mat_units[0]
    fs = unit_0.Meta_data.parameters.Sample_Rate
    trial_key = dict(session_info, trial_counts=len(unit_0.Trial_info.Trial_types))
    if trial_key not in acquisition.TrialSet.proj():
        print('\nInsert trial information')
        acquisition.TrialSet.insert1(trial_key, allow_direct_insert=True, ignore_extra_fields = True)

        for tr_idx in tqdm(np.arange(len(unit_0.Trial_info.Trial_types)) + 1):
            trial_key['trial_id'] = tr_idx
            trial_key['start_time'] = (unit_0.Trial_info.Trial_onset_vector_original[tr_idx] / fs
                                       if 'Trial_onset_vector_original' in unit_0.Trial_info._fieldnames else None)
            trial_key['stop_time'] = None  # hard-coded here, no trial-end times found in data
            trial_key['trial_stim_present'] = unit_0.Behavior.stim_trial_vector[tr_idx] != 0
            trial_key['trial_is_good'] = unit_0.Trial_info.Trial_range_to_analyze[0] <= tr_idx <= unit_0.Trial_info.Trial_range_to_analyze[-1]
            trial_key['trial_type'], trial_key['trial_response'] = trial_type_and_response_dict[
                unit_0.Behavior.Trial_types_of_response_vector[tr_idx]]
            trial_key['delay_duration'] = unit_0.Behavior.delay_dur[tr_idx]
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

    # ==================== behavioral and photostim ====================


# ====================== Starting import and compute procedure ======================
print('======== Populate() Routine =====')


