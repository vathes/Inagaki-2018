'''
Schema of behavioral information.
'''
import re
import os
from datetime import datetime
import sys

import numpy as np
import scipy.io as sio
import datajoint as dj
import h5py as h5

from . import utilities, acquisition, analysis, intracellular
from . import intracellular_path

schema = dj.schema(dj.config['custom'].get('database.prefix', '') + 'behavior')
sess_data_dir = os.path.join(intracellular_path, 'Data')


@schema
class TrialSegmentedLickTrace(dj.Imported):
    definition = """
    -> acquisition.TrialSet.Trial
    -> analysis.TrialSegmentationSetting
    ---
    segmented_lick_left_on: longblob  # (s), lick left onset times (based on contact of lick port)
    segmented_lick_left_off: longblob  # (s), lick left offset times (based on contact of lick port)
    segmented_lick_right_on: longblob  # (s), lick right onset times (based on contact of lick port)
    segmented_lick_right_off: longblob  # (s), lick right offset times (based on contact of lick port)
    """

    key_source = ((acquisition.TrialSet & (acquisition.Session.ExperimentType
                                           & {'experiment_type': 'intracellular'}))
                  * analysis.TrialSegmentationSetting)

    def make(self, key):
        # ============ Dataset ============
        # Get the Session definition from the keys of this session
        sess_data_file = utilities.find_session_matched_matfile(
            sess_data_dir, dict(key, cell_id=(intracellular.Cell & key).fetch1('cell_id')))
        if sess_data_file is None:
            raise FileNotFoundError(f'Intracellular import failed: ({key["subject_id"]} - {key["session_time"]})')
        mat_data = sio.loadmat(sess_data_file, struct_as_record = False, squeeze_me = True)['wholeCell']
        # get event, pre/post stim duration
        event_name, pre_stim_dur, post_stim_dur = (analysis.TrialSegmentationSetting & key).fetch1(
            'event', 'pre_stim_duration', 'post_stim_duration')
        pre_stim_dur = float(pre_stim_dur)
        post_stim_dur = float(post_stim_dur)
        #  ============= Now read the data and start ingesting =============
        self.insert({**key,
                     **get_single_trial_lick_times(k, mat_data, event_name, pre_stim_dur, post_stim_dur)}
                    for k in (acquisition.TrialSet.Trial & key).fetch('KEY'))
        print(f'Perform trial-segmentation of lick traces for trial: {key["session_id"]}')


def get_single_trial_lick_times(trial_key, mat_data, event_name, pre_stim_dur, post_stim_dur):
    lick_times = {k_n: getattr(mat_data.behavioral_data.behav_timing[trial_key['trial_id'] - 1], n)
                  for k_n, n in zip(['segmented_lick_left_on', 'segmented_lick_left_off',
                                     'segmented_lick_right_on', 'segmented_lick_right_off'],
                                    ['lickL_on_time', 'lickL_off_time',
                                     'lickR_on_time', 'lickR_off_time'])}
    # get event time
    try:
        event_time_point = analysis.get_event_time(event_name, trial_key)
    except analysis.EventChoiceError as e:
        raise e
    # check if pre/post stim dur is within start/stop time
    trial_start, trial_stop = (acquisition.TrialSet.Trial & trial_key).fetch1('start_time', 'stop_time')
    if trial_start and event_time_point - pre_stim_dur < 0:
        print('Warning: Out of bound pre-stim dur, select from start-time (t=0)')
        pre_stim_dur = event_time_point
    if trial_stop and event_time_point + post_stim_dur > trial_stop:
        print('Warning: Out of bound post-stim dur, set to trial end time')
        post_stim_dur = trial_stop - event_time_point

    for k, v in lick_times.items():
        v = np.array([v]) if isinstance(v, (float, int)) else v
        trial_key[k] = v[np.logical_and((v >= (event_time_point - pre_stim_dur)),
                                        (v <= (event_time_point + post_stim_dur)))] - event_time_point
    return trial_key

