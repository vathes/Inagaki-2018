from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)

# ====================== Starting import and compute procedure ======================
print('======== Populate() Intracellular Routine =====')
intracellular.MembranePotential.populate(reserve_jobs=True, suppress_errors=True)
intracellular.CurrentInjection.populate(reserve_jobs=True, suppress_errors=True)
intracellular.CellSpikeTimes.populate(reserve_jobs=True, suppress_errors=True)

intracellular.TrialSegmentedMembranePotential.populate(reserve_jobs=True, suppress_errors=True)
intracellular.TrialSegmentedCurrentInjection.populate(reserve_jobs=True, suppress_errors=True)
intracellular.TrialSegmentedCellSpikeTimes.populate(reserve_jobs=True, suppress_errors=True)

behavior.TrialSegmentedLickTrace.populate(reserve_jobs=True, suppress_errors=True)

print('======== Populate() Extracellular Routine =====')
analysis.RealignedEvent.populate(reserve_jobs=True, suppress_errors=True)
extracellular.UnitSpikeTimes.populate(reserve_jobs=True, suppress_errors=True)
extracellular.TrialSegmentedUnitSpikeTimes.populate(reserve_jobs=True, suppress_errors=True)

