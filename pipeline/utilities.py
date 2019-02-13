import os
from datetime import datetime
import re

import h5py as h5
import numpy as np

from . import reference, acquisition


# datetime format - should probably read this from a config file and not hard coded here
datetimeformat_ymdhms = '%Y%m%d'
datetimeformat_ymd = '%m%d%y'

def parse_prefix(line):
    cover = len(datetime.now().strftime(datetimeformat_ymdhms))
    try:
        return datetime.strptime(line[:cover], datetimeformat_ymdhms)
    except Exception as e:
        msg = f'Error:  {str(e)} \n'
        cover = len(datetime.now().strftime(datetimeformat_ymd))
        try:
            return datetime.strptime(line[:cover], datetimeformat_ymd)
        except Exception as e:
            print(f'{msg}\t{str(e)}\n\tReturn None')
            return None    


def find_session_matched_matfile(sess_data_dir, key):
        ############## Dataset #################
        sess_data_file = None
        which_data = re.search('Probe|WholeCell', sess_data_dir).group()

        # Search the filenames to find a match for "this" session (based on key)
        if which_data == 'WholeCell':
            cell_id = key['cell_id']
            exp_type = re.search('regular|EPSP',
                                 ';'.join((acquisition.Session.ExperimentType & key).fetch(
                                     'experiment_type'))).group()
            fnames = (f for f in os.listdir(sess_data_dir) if f.find('.mat') != -1)
            for f in fnames:
                if f.find(f'{cell_id}_{exp_type}') != -1:
                    sess_data_file = os.path.join(sess_data_dir, f)
                    break
        elif which_data == 'Probe':
            fnames = np.hstack([os.path.join(dir_files[0], f) for f in dir_files[2] if f.find('.mat') != -1]
                               for dir_files in os.walk(sess_data_dir) if len(dir_files[1]) == 0)
            for f in fnames:
                if f.find(key['session_id']) != -1:
                    sess_data_file = f
                    break

        # If session not found from dataset, break
        if sess_data_file is None:
            print(f'Session not found! - Subject: {key["subject_id"]} - Date: {key["session_time"]}')
            return None
        else: 
            print(f'Found datafile: {sess_data_file}')
            return sess_data_file
 
       
def get_brain_hemisphere(brain_region):
    # hemisphere: left-hemisphere is ipsi, so anything contra is right
    if re.search('Contra\s?', brain_region) is not None:
        hemi = 'right'
        brain_region = re.sub('Contra\s?', '', brain_region)
    else:
        hemi = 'left'
    return brain_region, hemi
        
        