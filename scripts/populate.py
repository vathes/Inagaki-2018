from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)

# prioritized_sessions = [{'session_id': 'Whole_cell_96_regular'},
#                         {'session_id': 'Whole_cell_96_EPSP'},
#                         {'session_id': 'HI127_031617'},
#                         {'session_id': 'HI152_060218'}]

prioritized_sessions = acquisition.Session()
settings = dict(reserve_jobs=True, suppress_errors=True)

# ====================== Starting import and compute procedure ======================
print('======== Populate() Intracellular Routine =====')
intracellular.MembranePotential.populate(prioritized_sessions, **settings)
intracellular.CurrentInjection.populate(prioritized_sessions, **settings)
intracellular.CellSpikeTimes.populate(prioritized_sessions, **settings)

intracellular.TrialSegmentedMembranePotential.populate(prioritized_sessions, **settings)
intracellular.TrialSegmentedCurrentInjection.populate(prioritized_sessions, **settings)
intracellular.TrialSegmentedCellSpikeTimes.populate(prioritized_sessions, **settings)

behavior.TrialSegmentedLickTrace.populate(prioritized_sessions, **settings)

print('======== Populate() Extracellular Routine =====')
analysis.RealignedEvent.populate(**settings)
extracellular.UnitSpikeTimes.populate(prioritized_sessions, **settings)
extracellular.TrialSegmentedUnitSpikeTimes.populate(prioritized_sessions, **settings)

