"""
qc_state_repair
===============

pDynamo's DFT/HF models (QCModelDFT) keep the fitting basis in
qcState.fitBases, an attribute that QCModelBaseState does not declare as
attributable, so it is not pickled. A DFT system that went through pickle
(a .pkl or .easy file, Clone(), a multiprocessing.Pool) raises

    AttributeError: 'QCModelBaseState' object has no attribute 'fitBases'

at the next energy. Setting fitBases alone is not enough (the C state of
the model is gone too: the process crashes), so the QC model is defined
again on the same QC atoms, and the NB model after it.

    repair_qc_state(system)   # returns True when the system was repaired
"""


def needs_repair(system):
    model = getattr(system, "qcModel", None)
    if model is None or not hasattr(model, "fitBasis"):
        return False
    state = getattr(system, "qcState", None)
    return state is not None and not hasattr(state, "fitBases")


def repair_qc_state(system):
    if not needs_repair(system):
        return False
    from pCore import Clone, Selection
    model = Clone(system.qcModel)
    qc_atoms = list(system.qcState.pureQCAtoms)
    nb_model = system.nbModel
    n_atoms = len(system.atoms)
    if len(qc_atoms) == n_atoms:
        system.DefineQCModel(model)
    else:
        system.DefineQCModel(model, qcSelection=Selection.FromIterable(qc_atoms))
    if nb_model is not None:
        system.DefineNBModel(Clone(nb_model))
    return True
