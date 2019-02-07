'''
Schema of session information.
'''
import datajoint as dj
from pipeline import acquisition

schema = dj.schema(dj.config.get('database.prefix', '') + 'ephys')

