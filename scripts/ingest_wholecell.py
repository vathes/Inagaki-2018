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
path = os.path.join('.', 'data', 'WholeCellData', 'Data')
fnames = os.listdir(path)

xlsname = 'SI_table_1_wc_cell_list.xlsx'
meta_data = pd.read_excel(os.path.join('.', 'data', 'WholeCellData', xlsname),
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
fnames = (f for f in os.listdir(path) if f.find('.mat') != -1)
for fname in fnames:
    mat_data = sio.loadmat(os.path.join(path, fname), struct_as_record = False, squeeze_me = True)['wholeCell']
    this_sess = meta_data.loc[f'Cell {mat_data.cell_id}']
    print(f'\nReading: {fname}')

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
                        session_id=uuid.uuid4().hex,
                        session_time=session_time)

    experimenters = ['Hidehiko Inagaki']  # hard-coded here
    experiment_types = this_sess.experiment_type
    experiment_types = [experiment_types] if isinstance(experiment_types, str) else experiment_types
    experiment_types.append(f'intracellular_{re.search("regular|EPSP", fname).group()}')

    # experimenter and experiment type (possible multiple experimenters or types)
    # no experimenter info
    acquisition.ExperimentType.insert(zip(experiment_types), skip_duplicates=True)

    if session_info not in acquisition.Session.proj():
        with acquisition.Session.connection.transaction:
            acquisition.Session.insert1(session_info, ignore_extra_fields=True)
            acquisition.Session.Experimenter.insert((dict(session_info, experiment_type=k) for k in experimenters), ignore_extra_fields=True)
            acquisition.Session.ExperimentType.insert((dict(session_info, experiment_type=k) for k in experiment_types), ignore_extra_fields=True)
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
    if brain_location not in reference.BrainLocation.proj():
        reference.BrainLocation.insert1(brain_location)

    # -- Cell
    cell_id = f'cell_{mat_data.cell_id}'
    cell_key = dict({**session_info, **brain_location},
                    cell_id=cell_id,
                    cell_type=mat_data.Pyr_or_GABA,
                    cell_depth=float(this_sess.depth_um),
                    device_name=ie_device)

    if cell_key not in acquisition.Cell.proj():
        acquisition.Cell.insert1(cell_key, ignore_extra_fields=True)

    # ==================== Trials ====================
    fs = mat_data.recording_data.sample_rate
    trial_key = dict(session_info, trial_counts=len(mat_data.behavioral_data.behav_timing))
    if trial_key not in acquisition.TrialSet.proj():
        print('\nInsert trial information')
        acquisition.TrialSet.insert1(trial_key, allow_direct_insert=True, ignore_extra_fields = True)

        for tr_idx, events_time in tqdm(enumerate(mat_data.behavioral_data.behav_timing)):
            trial_key['trial_id'] = tr_idx
            trial_key['start_time'] = mat_data.behavioral_data.trial_onset_bin[tr_idx] / fs
            trial_key['stop_time'] = trial_key['start_time'] + events_time.end_time
            trial_key['trial_stim_present'] = bool(mat_data.behavioral_data.AOM_on_or_off[tr_idx])
            trial_key['trial_is_good'] = True  #  no info of trial good/bad status, assuming all trials are good
            trial_key['trial_type'], trial_key['trial_response'] = trial_type_and_response_dict[
                mat_data.behavioral_data.trial_type_vector[tr_idx]]
            trial_key['delay_duration'] = 1.2  # hard-coded here (the same for whole cell)
            acquisition.TrialSet.Trial.insert1(trial_key, ignore_extra_fields = True, skip_duplicates = True,
                                               allow_direct_insert = True)

            # ======== Now add trial event timing to the EventTime part table ====
            # -- events timing
            acquisition.TrialSet.EventTime.insert((dict(trial_key, trial_event=k, event_time=dict(
                events_time.__dict__,
                trial_start=trial_key['start_time'],
                trial_stop=trial_key['stop_time'],
                first_lick = min([events_time.__getattribute__(l)[0]
                                  for l in ('lickL_on_time', 'lickR_on_time')
                                  if events_time.__getattribute__(l).size > 0]),
                current_injection_start=mat_data.behavioral_data.tail_current_injection_onset_bin[tr_idx] / fs)[k])
                                                   for k in ['trial_start', 'trial_stop', 'cue_start',
                                                             'cue_end', 'sampling_start', 'delay_start',
                                                             'current_injection_start']),
                                                  ignore_extra_fields = True, skip_duplicates = True,
                                                  allow_direct_insert = True)

    # ==================== behavioral and photostim ====================


# ====================== Starting import and compute procedure ======================
print('======== Populate() Routine =====')
# -- Intracellular
acquisition.IntracellularAcquisition.populate()

