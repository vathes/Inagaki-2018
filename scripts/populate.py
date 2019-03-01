from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)

prioritized_sessions = [{'session_id': 'Whole_cell_96_regular'},
                        {'session_id': 'Whole_cell_96_EPSP'},
                        {'session_id': 'HI127_031617'},
                        {'session_id': 'HI152_060218'}]

# ====================== Starting import and compute procedure ======================
print('======== Populate() Intracellular Routine =====')
intracellular.MembranePotential.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)
intracellular.CurrentInjection.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)
intracellular.CellSpikeTimes.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)

intracellular.TrialSegmentedMembranePotential.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)
intracellular.TrialSegmentedCurrentInjection.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)
intracellular.TrialSegmentedCellSpikeTimes.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)

behavior.TrialSegmentedLickTrace.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)

print('======== Populate() Extracellular Routine =====')
analysis.RealignedEvent.populate(reserve_jobs=True, suppress_errors=True)
extracellular.UnitSpikeTimes.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)
extracellular.TrialSegmentedUnitSpikeTimes.populate(prioritized_sessions, reserve_jobs=True, suppress_errors=True)

