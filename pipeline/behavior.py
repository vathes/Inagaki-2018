'''
Schema of session information.
'''
import re
import os
from datetime import datetime

import numpy as np
import scipy.io as sio
import datajoint as dj

from . import acquisition, reference
from .utilities import get_one_from_nested_array, get_list_from_nested_array, datetimeformat_ydm, datetimeformat_ymd

schema = dj.schema(dj.config.get('database.prefix', '') + 'behavior')
       
