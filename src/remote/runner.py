#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Runs one EasyHybrid job on the remote machine:
#      python3 -m remote.runner <job folder>
#  (PYTHONPATH holds the copy of src/ and pDynamo.) The job folder has
#  system.pkl, parameters.pkl and job.json; the results go to output/,
#  exactly as a local run writes them to its folder, and status.json says
#  how it ended: {"state": "finished" | "error", "logfile": <relative to
#  output/>, "warning", "error", "traceback", "final_coordinates"}.
#  No GTK: only pDynamo and src/pdynamo.

import json
import os
import pickle
import socket
import sys
import time
import traceback


def _write_status(job_dir, **data):
    data.setdefault("time", time.strftime("%Y-%m-%d %H:%M:%S"))
    data.setdefault("host", socket.gethostname())
    tmp = os.path.join(job_dir, "status.json.tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, indent=1, default=str)
    os.replace(tmp, os.path.join(job_dir, "status.json"))


def _simulation_classes():
    from pdynamo import p_methods as pm
    return {
        "Energy_Single_Point": pm.EnergyCalculation,
        "Energy_Refinement": pm.EnergyRefinement,
        "Geometry_Optimization": pm.GeometryOptimization,
        "Molecular_Dynamics": pm.MolecularDynamics,
        "Relaxed_Surface_Scan": pm.RelaxedSurfaceScan,
        "Advanced_Relaxed_Surface_Scan": pm.AdvancedRelaxedSurfaceScan,
        "Umbrella_Sampling": pm.UmbrellaSampling,
        "Nudged_Elastic_Band": pm.ChainOfStatesOptimizePath,
        "Normal_Modes": pm.NormalModes,
        "Conjugate_Peak_Refinement": pm.ConjugatePeakRefinementRun,
    }


def _logfile_guess(parameters, job_id):
    """ Same first guess as simulations_mixin._configure_logfile (the runners
        correct it to their real path). """
    folder = parameters["folder"]
    if parameters.get("filename"):
        return os.path.join(folder, parameters["filename"] + ".log")
    if parameters.get("trajectory_name"):
        return os.path.join(folder, parameters["trajectory_name"], "output.log")
    if parameters.get("traj_folder_name"):
        return os.path.join(folder, parameters["traj_folder_name"] + ".ptGeo", "output.log")
    return os.path.join(folder, "{}.log".format(job_id))


_ENGINES = {"QCModelORCA": ("ORCA", "PDYNAMO3_ORCACOMMAND", "orca"),
            "QCModelXTB": ("xTB", "PDYNAMO3_XTBCOMMAND", "xtb"),
            "QCModelDFTB": ("DFTB+", "PDYNAMO3_DFTBCOMMAND", "dftb")}


def qc_engine(system):
    """ (label, variable, profile key) of the external QC program of a system, or None. """
    model = getattr(system, "qcModel", None)
    return _ENGINES.get(type(model).__name__) if model is not None else None


def _prepare_qc(system, job_dir):
    """ External QC programs (ORCA, xTB, DFTB+) keep the scratch paths of the
        computer the system came from, in the model and in its state: they
        move to a folder of this job (EASYHYBRID_QC_SCRATCH, e.g. a node-local
        disk, overrides it), which also keeps simultaneous jobs apart (the
        parallel runners add process_<i> under it). DFTB+ parameters sent
        with the job (input/skf) replace the local skfPath. The program must
        be found here: a clear error otherwise. """
    engine = qc_engine(system)
    if engine is None:
        return None
    label, variable, key = engine
    model, state = system.qcModel, getattr(system, "qcState", None)
    base = os.environ.get("EASYHYBRID_QC_SCRATCH") or job_dir
    scratch = os.path.join(base, "qc_scratch")
    os.makedirs(scratch, exist_ok=True)
    model.scratch = scratch
    if state is not None and hasattr(state, "DeterminePaths"):
        state.DeterminePaths(scratch, deleteJobFiles=getattr(state, "deleteJobFiles", True))
    shipped = os.path.join(job_dir, "input", "skf")
    if hasattr(model, "skfPath") and os.path.isdir(shipped):
        model.skfPath = shipped
    try:
        model.command
    except Exception:
        raise RuntimeError("{} was not found on {}: give its path in the remote machine's profile "
                           "(QC programs) or set {} in the setup lines. ({}={})".format(
                               label, socket.gethostname(), variable, variable, os.environ.get(variable)))
    if label == "ORCA":
        import re
        import shutil
        text = " ".join(str(k) for k in (getattr(model, "keywords", None) or []))
        if (re.search(r"\bPAL(\d+)\b", text, re.I) or re.search(r"nprocs\s+\d+", text, re.I)) and not shutil.which("mpirun"):
            raise RuntimeError("Parallel ORCA (PAL) needs mpirun, not found on {}: load the MPI module ORCA was "
                               "built with in the setup lines (e.g. module load openmpi/4.1.x).".format(socket.gethostname()))
    print("QC program: {} = {}, scratch {}".format(label, model.command, scratch))
    return label


def _qc_error_hint(job_dir):
    """ The end of the QC program's own log when it failed (pDynamo renames
        the files to error_*): the real reason, rather than the Python error
        that follows it. """
    import glob
    logs = sorted(glob.glob(os.path.join(job_dir, "qc_scratch", "**", "error_*.log"), recursive=True) +
                  glob.glob(os.path.join(os.environ.get("EASYHYBRID_QC_SCRATCH") or job_dir, "qc_scratch", "**",
                                         "error_*.out"), recursive=True))
    if not logs:
        return ""
    keys = ("error", "not enough slots", "command not found", "cannot", "abort", "fail", "not found", "illegal")
    lines = [l.rstrip() for l in open(logs[0], errors="replace") if l.strip()]
    found = [l for l in lines if any(k in l.lower() for k in keys)]
    return "\n--- {} ---\n{}".format(os.path.relpath(logs[0], job_dir), "\n".join((found or lines)[-12:]))


def run(job_dir):
    job_dir = os.path.abspath(job_dir)
    with open(os.path.join(job_dir, "job.json")) as f:
        job = json.load(f)
    _write_status(job_dir, state="running", started=time.strftime("%Y-%m-%d %H:%M:%S"))
    output = os.path.join(job_dir, "output")
    os.makedirs(output, exist_ok=True)
    try:
        from pBabel import ImportSystem
        from pCore import Pickle
        with open(os.path.join(job_dir, "parameters.pkl"), "rb") as f:
            parameters = pickle.load(f)
        system = ImportSystem(os.path.join(job_dir, "system.pkl"))
        _prepare_qc(system, job_dir)
        try:
            from util.qc_state_repair import repair_qc_state
            repair_qc_state(system)
        except Exception as error:
            print("Could not check the QC state of the system:", error)
        parameters["system"] = system
        parameters["folder"] = output
        if parameters.get("source_folder"):
            parameters["source_folder"] = os.path.join(job_dir, parameters["source_folder"])
        parameters["queue"] = None
        parameters["logfile"] = _logfile_guess(parameters, job["id"])
        sim_type = parameters["simulation_type"]
        cls = _simulation_classes().get(sim_type)
        if cls is None:
            raise ValueError("Simulation type {} cannot run remotely.".format(sim_type))
        if sim_type == "Energy_Single_Point":
            parameters["energy"] = True
        cls().run(parameters)
        final = os.path.join(output, "final_coordinates.pkl")
        Pickle(final, system.coordinates3)
        logfile = parameters.get("logfile")
        _write_status(job_dir, state="finished", started=job.get("started"),
                      logfile=os.path.relpath(logfile, output) if logfile else None,
                      warning=parameters.get("warning"), final_coordinates="final_coordinates.pkl")
        return 0
    except BaseException as error:
        _write_status(job_dir, state="error", error="{}: {}{}".format(type(error).__name__, error, _qc_error_hint(job_dir)),
                      traceback=traceback.format_exc())
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(run(sys.argv[1]))
