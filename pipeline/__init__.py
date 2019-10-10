import os
import datajoint as dj
import pathlib


if 'custom' not in dj.config:
    raise KeyError('"custom" portion of the dj_local_conf.json not found, see README for instruction')

intracellular_path = pathlib.Path(dj.config['custom'].get('intracellular_directory')).as_posix()
extracellular_path = pathlib.Path(dj.config['custom'].get('extracellular_directory')).as_posix()
