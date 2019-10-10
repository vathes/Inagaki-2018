from pipeline import (reference, subject, acquisition, stimulation, analysis,
                      intracellular, extracellular, behavior, utilities)

settings = dict(reserve_jobs=True, suppress_errors=True)

# ====================== Starting import and compute procedure ======================
print('======== Populate() Intracellular Routine =====')
intracellular.MembranePotential.populate(**settings)
intracellular.CurrentInjection.populate(**settings)
intracellular.CellSpikeTimes.populate(**settings)

intracellular.TrialSegmentedMembranePotential.populate(**settings)
intracellular.TrialSegmentedCurrentInjection.populate(**settings)
intracellular.TrialSegmentedCellSpikeTimes.populate(**settings)

behavior.TrialSegmentedLickTrace.populate(**settings)

print('======== Populate() Extracellular Routine =====')
analysis.RealignedEvent.populate(**settings)
extracellular.UnitSpikeTimes.populate(**settings)
extracellular.TrialSegmentedUnitSpikeTimes.populate(**settings)
