"""
Testes dos jobs remotos (src/remote): perfis, script do SLURM, leitura dos
estados do escalonador, e o ciclo completo enviar -> rodar -> acompanhar ->
buscar com o transporte local (sem ssh) e o escalonador "direct". Tudo em
pastas temporarias: nada em ~/.easyhybrid.
"""
import os
import time

import pytest

from remote import profiles, schedulers

PDYNAMO_ENV = "/home/fernando/programs/pDynamo3/installation/shellScripts/environment_bash.com"


def test_profiles_roundtrip(tmp_path):
    path = str(tmp_path / "remotes.json")
    with pytest.raises(ValueError):
        profiles.save({"name": "x"}, path)                      # . no host, no folder
    saved = profiles.save({"name": "hpc", "host": "hpc", "workdir": "/scratch/u/eh", "scheduler": "slurm",
                           "slurm": {"partition": "cpu"}}, path)
    assert saved["slurm"]["time"] == "24:00:00" and saved["slurm"]["partition"] == "cpu"
    assert profiles.get("hpc", path)["workdir"] == "/scratch/u/eh"
    assert profiles.validate(dict(saved, scheduler="pbs"))
    assert profiles.remove("hpc", path) and profiles.load_all(path) == {}


def test_slurm_script():
    profile = profiles._merge({"name": "hpc", "host": "hpc", "workdir": "/s", "scheduler": "slurm",
                               "setup": "module load python/3.10\nsource /opt/pDynamo3/env.sh",
                               "slurm": {"partition": "computeq", "account": "lab", "mem": "4G",
                                         "extra": "--qos=normal"}})
    text = schedulers.get("slurm").script(profile, {"id": "j1", "name": "my scan", "cpus": 8},
                                          "/s/jobs/j1", "/s/easyhybrid/src")
    for line in ("#!/bin/bash -l", "#SBATCH --cpus-per-task=8", "#SBATCH --partition=computeq",
                 "#SBATCH --account=lab", "#SBATCH --mem=4G", "#SBATCH --qos=normal", "#SBATCH --job-name=eh_my_scan",
                 "module load python/3.10", "export PYTHONPATH=/s/easyhybrid/src", "-m remote.runner /s/jobs/j1"):
        assert line in text, line


class _FakeTransport:
    def __init__(self, answers):
        self.answers, self.commands = answers, []

    def run(self, command, timeout=120, check=True, stdin=None):
        self.commands.append(command)
        for key, out in self.answers.items():
            if command.startswith(key):
                return 0, out, ""
        return 0, "", ""


def test_slurm_states():
    transport = _FakeTransport({"squeue": "101 RUNNING\n102 PENDING\n", "sacct": "103|COMPLETED\n104|CANCELLED by 5\n"})
    states = schedulers.get("slurm").states(transport, ["101", "102", "103", "104", "105"])
    assert states == {"101": "RUNNING", "102": "PENDING", "103": "COMPLETED", "104": "CANCELLED", "105": "GONE"}
    assert "sacct" in transport.commands[1] and "103,104,105" in transport.commands[1]


def test_remote_cycle_local(tmp_path, monkeypatch):
    """ A 1D scan of bALA sent to "this computer" through the remote path. """
    pBabel = pytest.importorskip("pBabel")
    if not os.path.exists(PDYNAMO_ENV):
        pytest.skip("pDynamo environment script not found")
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelFull
    from remote import jobs
    monkeypatch.setattr(jobs, "CONTROL_DIR", str(tmp_path / "control"))
    pf = str(tmp_path / "remotes.json")
    profiles.save({"name": "here", "host": "localhost-test", "workdir": str(tmp_path / "remote_root"),
                   "setup": "source " + PDYNAMO_ENV, "scheduler": "direct"}, pf)
    manager = jobs.RemoteJobs(registry=str(tmp_path / "jobs.json"), profiles_file=pf)
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    system = pBabel.ImportSystem(os.path.join(repo, "examples", "tutorials", "01_bala", "bALA_not_opt.mol"))
    system.DefineMMModel(MMModelOPLS.WithParameterSet("protein"))
    system.DefineNBModel(NBModelFull.WithDefaults())
    system.e_id, system.e_tag, system.e_step_counter = 0, "bALA", 0
    system.e_treeview_iter = object()                     # . not picklable: must be left out of system.pkl
    params = {"simulation_type": "Relaxed_Surface_Scan", "logFrequency": 50, "optimizer": "ConjugatedGradient",
              "maximumIterations": 1000.0, "rmsGradientTolerance": 0.1, "traj_folder_name": "scan_r",
              "vobject_name": "scan_r", "_is_ts_centered": False, "RC2": None, "NmaxThreads": 1,
              "traj_type": "pklfolder",
              "RC1": {"rc_type": "simple_distance", "ATOMS": [5, 17], "ATOM_NAMES": ["O5", "H17"], "dminimum": 2.1,
                      "sigma_pk1pk3": None, "sigma_pk3pk1": None, "nsteps": 4, "dincre": 0.2, "force_constant": 4000.0}}
    local = tmp_path / "results"
    local.mkdir()
    job = manager.prepare(system, params, "here", str(local))
    assert job["state"] == "preparing" and isinstance(system.e_treeview_iter, object)
    job = manager.send(job["id"])
    assert job["state"] == "running" and job["scheduler_id"].isdigit()
    for _ in range(90):
        time.sleep(1)
        job = manager.refresh([job["id"]])[0]
        if job["state"] in jobs.FINAL_STATES:
            break
    assert job["state"] == "finished", job
    job = manager.fetch(job["id"])
    assert job["fetched"] and job["logfile"] == str(local / "scan_r.ptGeo" / "output.log")
    assert os.path.isfile(job["logfile"]) and os.path.isfile(job["final_coordinates"])
    assert len([f for f in os.listdir(local / "scan_r.ptGeo") if f.endswith(".pkl")]) == 4
    assert not os.path.exists(tmp_path / "control" / "bundles" / job["id"])        # . bundle removed


def test_runner_reports_errors(tmp_path):
    """ A job that fails writes state "error" with the message. """
    import json
    from remote import runner
    (tmp_path / "job.json").write_text(json.dumps({"id": "bad"}))
    assert runner.run(str(tmp_path)) == 1
    status = json.loads((tmp_path / "status.json").read_text())
    assert status["state"] == "error" and "parameters.pkl" in status["error"]


def test_session_remote_flow(tmp_path, monkeypatch):
    """ pRemote (the session side): run_simulation() with a remote target ->
        Process Manager row -> polling -> automatic fetch -> final object -> Finished. """
    pBabel = pytest.importorskip("pBabel")
    if not os.path.exists(PDYNAMO_ENV):
        pytest.skip("pDynamo environment script not found")
    gi = pytest.importorskip("gi")
    gi.require_version("Gtk", "3.0")
    from gi.repository import GLib
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelFull
    from remote import jobs
    from pdynamo.pDynamo2EasyHybrid import remote_mixin
    from pdynamo.pDynamo2EasyHybrid.remote_mixin import pRemote
    monkeypatch.setattr(jobs, "CONTROL_DIR", str(tmp_path / "control"))
    monkeypatch.setattr(jobs, "REGISTRY_FILE", str(tmp_path / "jobs.json"))
    monkeypatch.setattr(profiles, "PROFILES_FILE", str(tmp_path / "remotes.json"))
    monkeypatch.setattr(remote_mixin, "POLL_SECONDS", 1)
    profiles.save({"name": "here", "host": "localhost-test", "workdir": str(tmp_path / "remote_root"),
                   "setup": "source " + PDYNAMO_ENV, "scheduler": "direct"})

    calls = []

    class Stub:
        def __getattr__(self, name):
            return lambda *a, **k: calls.append((name, a, k)) or 0

    class Main:
        process_manager_window = Stub()
        bottom_notebook = Stub()
        simple_dialog = Stub()
        job_history_liststore = {0: [None] * 10}

    class Session(pRemote):
        def __init__(self, system):
            self.main, self.psystem, self.active_id = Main(), {0: system}, 0
            self.vm_session = type("V", (), {"gen_random_tag_string": lambda self, length=6: "abcdef"})()
            self.objects = []

        def get_hamiltonian_type(self, e_id):
            return "MM"

        def _apply_restraints(self, parameters):
            return parameters

        def _add_vismol_object_to_easyhybrid_session(self, system, name, key6):
            self.objects.append((name, key6, len(system.atoms)))
            return type("O", (), {})()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    system = pBabel.ImportSystem(os.path.join(repo, "examples", "tutorials", "01_bala", "bALA_not_opt.mol"))
    system.DefineMMModel(MMModelOPLS.WithParameterSet("protein"))
    system.DefineNBModel(NBModelFull.WithDefaults())
    system.e_id, system.e_tag, system.e_step_counter, system.e_job_history = 0, "bALA", 0, {}
    session = Session(system)
    session.main.process_manager_window.add_new_process = lambda **kw: 0
    out = tmp_path / "results"
    out.mkdir()
    params = {"simulation_type": "Geometry_Optimization", "trajectory_name": None, "dialog": False,
              "folder": str(out), "optimizer": "ConjugatedGradient", "maximumIterations": 600, "logFrequency": 10,
              "save_frequency": 10, "rmsGradientTolerance": 0.1, "save_format": None, "save_traj": False,
              "save_pdb": False, "vobject_name": "opt", "remote": "here"}
    job = session.run_remote(params, "here")
    assert job and system.e_step_counter == 1 and system.e_job_history[0]["remote_job"] == job["id"]
    loop, deadline = GLib.MainLoop(), time.time() + 90

    def check():
        if system.e_job_history[0].get("status") in ("Finished", "Error!") or time.time() > deadline:
            loop.quit()
            return False
        session.remote_refresh_now([job["id"]])
        return True

    GLib.timeout_add(1000, check)
    loop.run()
    entry = system.e_job_history[0]
    assert entry["status"] == "Finished", (entry, calls)
    assert entry["logfile"] and os.path.isfile(entry["logfile"])
    assert session.objects == [("Geometry_Optimization_0", "abcdef", 22)]
    assert ("set_status", (0, "Finished"), {"warning": None}) in calls


def test_remote_reattach_after_restart(tmp_path, monkeypatch):
    """ pRemote.remote_reattach(): closing EasyHybrid while a remote job is
        running, then reopening the same session. build_liststore_from_job_history()
        rebuilds a static Process Manager row from the pickled system.e_job_history,
        but a brand-new session never submitted that job, so its _remote_rows is
        empty -- without reattach, _remote_show()/_remote_done() would see
        row is None and silently do nothing. remote_reattach() wires the
        reconstructed row back to the registry, so polling and the automatic
        fetch work again. """
    pBabel = pytest.importorskip("pBabel")
    if not os.path.exists(PDYNAMO_ENV):
        pytest.skip("pDynamo environment script not found")
    gi = pytest.importorskip("gi")
    gi.require_version("Gtk", "3.0")
    from gi.repository import GLib
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelFull
    from remote import jobs
    from pdynamo.pDynamo2EasyHybrid import remote_mixin
    from pdynamo.pDynamo2EasyHybrid.remote_mixin import pRemote
    monkeypatch.setattr(jobs, "CONTROL_DIR", str(tmp_path / "control"))
    monkeypatch.setattr(jobs, "REGISTRY_FILE", str(tmp_path / "jobs.json"))
    monkeypatch.setattr(profiles, "PROFILES_FILE", str(tmp_path / "remotes.json"))
    monkeypatch.setattr(remote_mixin, "POLL_SECONDS", 1)
    profiles.save({"name": "here", "host": "localhost-test", "workdir": str(tmp_path / "remote_root"),
                   "setup": "source " + PDYNAMO_ENV, "scheduler": "direct"})

    class Stub:
        def __getattr__(self, name):
            return lambda *a, **k: 0

    class Main:
        process_manager_window = Stub()
        bottom_notebook = Stub()
        simple_dialog = Stub()
        job_history_liststore = {0: [None] * 10}

    class Session(pRemote):
        def __init__(self, system):
            self.main, self.psystem, self.active_id = Main(), {0: system}, 0
            self.vm_session = type("V", (), {"gen_random_tag_string": lambda self, length=6: "abcdef"})()
            self.objects = []

        def get_hamiltonian_type(self, e_id):
            return "MM"

        def _apply_restraints(self, parameters):
            return parameters

        def _add_vismol_object_to_easyhybrid_session(self, system, name, key6):
            self.objects.append((name, key6, len(system.atoms)))
            return type("O", (), {})()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    system = pBabel.ImportSystem(os.path.join(repo, "examples", "tutorials", "01_bala", "bALA_not_opt.mol"))
    system.DefineMMModel(MMModelOPLS.WithParameterSet("protein"))
    system.DefineNBModel(NBModelFull.WithDefaults())
    system.e_id, system.e_tag, system.e_step_counter, system.e_job_history = 0, "bALA", 0, {}

    # . "session 1": the session that submits the job (run_simulation() with
    #   a remote target), like a user pressing Run before closing EasyHybrid.
    session1 = Session(system)
    out = tmp_path / "results"
    out.mkdir()
    params = {"simulation_type": "Geometry_Optimization", "trajectory_name": None, "dialog": False,
              "folder": str(out), "optimizer": "ConjugatedGradient", "maximumIterations": 600, "logFrequency": 10,
              "save_frequency": 10, "rmsGradientTolerance": 0.1, "save_format": None, "save_traj": False,
              "save_pdb": False, "vobject_name": "opt", "remote": "here"}
    job = session1.run_remote(params, "here")
    assert job and system.e_job_history[0]["remote_job"] == job["id"]
    time.sleep(1)               # . lets the submit thread's manager.send() actually run
    job_id = job["id"]

    # . "EasyHybrid is closed": the whole process (and every GLib source it
    #   registered) would be gone. This test runs in one process for both
    #   "sessions", so session1's own poll source -- on the SAME default
    #   GLib main context session2's loop below also runs on -- must be
    #   removed explicitly, or it would keep polling/fetching through
    #   session1 (which still has job_id in ITS _remote_rows) once that
    #   loop starts, racing session2's own fetch and making this test
    #   depend on which one wins instead of on remote_reattach().
    if getattr(session1, "_remote_poll_source", None) is not None:
        GLib.source_remove(session1._remote_poll_source)

    # . "session 2": a brand-new pDynamoSession after reopening the same .easy
    #   -- same system (e_job_history is pickled with it), but a fresh pRemote
    #   mixin state: no _remote_jobs_manager / _remote_rows yet.
    session2 = Session(system)
    assert session2.remote_reattach(0, 0, job_id, 0)
    assert job_id in session2._remote_rows

    loop, deadline = GLib.MainLoop(), time.time() + 90

    def check():
        if system.e_job_history[0].get("status") in ("Finished", "Error!") or time.time() > deadline:
            loop.quit()
            return False
        session2.remote_refresh_now([job_id])
        return True

    GLib.timeout_add(1000, check)
    loop.run()
    entry = system.e_job_history[0]
    assert entry["status"] == "Finished", entry
    assert entry["logfile"] and os.path.isfile(entry["logfile"])
    assert session2.objects == [("Geometry_Optimization_0", "abcdef", 22)]


def test_remote_reattach_unknown_job(tmp_path, monkeypatch):
    """ A row whose job is no longer in the registry (cleaned up, or the
        submission never actually reached it) is left as a static row:
        no crash, no polling started. """
    from remote import jobs
    from pdynamo.pDynamo2EasyHybrid.remote_mixin import pRemote
    monkeypatch.setattr(jobs, "CONTROL_DIR", str(tmp_path / "control"))
    monkeypatch.setattr(jobs, "REGISTRY_FILE", str(tmp_path / "jobs.json"))

    class Session(pRemote):
        def __init__(self):
            self.main = None

    session = Session()
    assert session.remote_reattach(0, 0, "no-such-job", 0) is False
    assert "no-such-job" not in session._remote_rows
    assert getattr(session, "_remote_poll_source", None) is None


def test_qc_requirements_and_orca_mpi_script():
    from remote.jobs import qc_requirements

    def system_with(name, **attrs):
        model = type(name, (), attrs)()
        return type("S", (), {"qcModel": model})()

    xtb = qc_requirements(system_with("QCModelXTB", parallel=4), {"NmaxThreads": 3})
    assert xtb["engine"] == "xTB" and xtb["cpus"] == 12
    orca = qc_requirements(system_with("QCModelORCA", keywords=["B3LYP 6-31G* PAL8"]), {})
    assert orca["engine"] == "ORCA" and orca["threads"] == 8 and orca["cpus"] == 8
    dftb = qc_requirements(system_with("QCModelDFTB", skfPath="/nonexistent"), {})
    assert dftb["engine"] == "DFTB+" and dftb["skf"] is None
    assert qc_requirements(type("S", (), {})(), {"NmaxThreads": 2})["engine"] is None
    profile = profiles._merge({"name": "h", "host": "h", "workdir": "/s", "scheduler": "slurm",
                               "qc_programs": {"orca": "/opt/orca/orca"}})
    text = schedulers.get("slurm").script(profile, {"id": "j", "cpus": 8, "qc_engine": "ORCA", "qc_threads": 8,
                                                    "processes": 1}, "/s/jobs/j", "/s/code")
    assert "#SBATCH --ntasks=8" in text and "#SBATCH --cpus-per-task=1" in text
    assert "export PDYNAMO3_ORCACOMMAND=/opt/orca/orca" in text and "binding_policy" not in text
    text = schedulers.get("slurm").script(profile, {"id": "j", "cpus": 8, "qc_engine": "ORCA", "qc_threads": 2,
                                                    "processes": 4}, "/s/jobs/j", "/s/code")
    assert "OMPI_MCA_hwloc_base_binding_policy=none" in text
    text = schedulers.get("slurm").script(profile, {"id": "j", "cpus": 4, "qc_engine": "xTB", "qc_threads": 4},
                                          "/s/jobs/j", "/s/code")
    assert "#SBATCH --ntasks=1" in text and "#SBATCH --cpus-per-task=4" in text


def test_runner_moves_qc_scratch(tmp_path, monkeypatch):
    from remote import runner

    class State:
        def DeterminePaths(self, scratch, deleteJobFiles=True):
            self.paths = {"Glob": os.path.join(scratch, "job")}

    class QCModelDFTB:
        scratch, skfPath = "/home/other/computer/dftbScratch", "/home/other/skf"

        @property
        def command(self):
            return "/bin/true"

    system = type("S", (), {"qcModel": QCModelDFTB(), "qcState": State()})()
    (tmp_path / "input" / "skf").mkdir(parents=True)
    assert runner._prepare_qc(system, str(tmp_path)) == "DFTB+"
    assert system.qcModel.scratch == str(tmp_path / "qc_scratch")
    assert system.qcState.paths["Glob"].startswith(str(tmp_path / "qc_scratch"))
    assert system.qcModel.skfPath == str(tmp_path / "input" / "skf")

    class QCModelXTB:
        scratch = "/elsewhere"

        @property
        def command(self):
            raise RuntimeError("not installed")

    monkeypatch.delenv("PDYNAMO3_XTBCOMMAND", raising=False)
    with pytest.raises(RuntimeError, match="xTB was not found"):
        runner._prepare_qc(type("S", (), {"qcModel": QCModelXTB(), "qcState": State()})(), str(tmp_path))


def test_runner_error_hint(tmp_path):
    from remote import runner
    (tmp_path / "qc_scratch").mkdir()
    (tmp_path / "qc_scratch" / "error_orcaJob.log").write_text(
        "banner\\nThere are not enough slots available in the system to satisfy the 2\\nmore\\n")
    hint = runner._qc_error_hint(str(tmp_path))
    assert "not enough slots" in hint and "error_orcaJob.log" in hint


def test_slurm_partitions():
    out = ("bigmem|up|28-00:00:00|4|202/566/0/768|(null)\n"
           "compute*|up|28-00:00:00|18|1654/1802/0/3456|(null)\n"
           "gpu|up|17-00:00:00|4|304/464/0/768|gpu:2,mps:200\n")
    rows = schedulers.get("slurm").partitions(_FakeTransport({"sinfo": out}))
    assert [r["name"] for r in rows] == ["compute", "bigmem", "gpu"]          # . default first
    assert rows[0]["default"] and rows[0]["cpus_idle"] == 1802 and rows[0]["cpus_total"] == 3456
    assert rows[2]["gres"] == "gpu:2,mps:200" and rows[1]["gres"] is None and rows[2]["time_limit"] == "17-00:00:00"


def test_discover_with_canned_answers(monkeypatch):
    """ Detect Settings on a cluster1-like machine: broken ORCA module, a module
        binary that needs its libraries, the Python that imports pDynamo. """
    from remote import discover as disc
    inventory = "\n".join([
        "HOST\tlog02", "USER\tu", "HOME\t/home/u", "NPROC\t32", "SYSPY\tPython 3.6.8",
        "CMD\tsbatch\t/usr/bin/sbatch", "CMD\tsinfo\t/usr/bin/sinfo", "HAS_MODULE\t1",
        "SCRATCH\t/scratch/u\t900000", "SCRATCH\t/home/u\t1000",
        "MODULE\tpython/3.10.13/gcc.8.5.0", "MODULE\tpython/3.12.1/gcc.8.5.0", "MODULE\torca/6.1.1",
        "MODULE\txtb/6.7.0", "MODULE\topenmpi/4.0.7/gcc.8.5.0/noncuda", "MODULE\topenmpi/4.1.8/gcc.8.5.0/noncuda",
        "MODPATH\torca/6.1.1\t/opt/apps/orca/6.1.1/bin", "MODPATH\txtb/6.7.0\t/opt/apps/xtb/6.7.0/bin",
        "PDYNAMO\t/home/u/programs/pDynamo3",
        "PROG\t/home/u/programs/xtb-6.7.1/bin/xtb",
        "PART\tcompute*|up|28-00:00:00|18|1654/1802/0/3456|(null)", "PART\tcomputeq|up|28-00:00:00|60|1/2/0/3|(null)",
        "ACCOUNT\tlab", "DEFACCOUNT\tlab"]) + "\n"

    class Fake:
        calls = 0

        def run(self, command, timeout=120, check=True, stdin=None):
            Fake.calls += 1
            if Fake.calls == 1:
                return 0, inventory, ""
            assert "python/3.10.13" in stdin and "/opt/apps/orca/6.1.1/orca" in stdin
            return 0, "\n".join([
                "COMMIT\t/home/u/programs/pDynamo3\tabc123",
                "ENVSCRIPT\t/home/u/programs/pDynamo3\t/home/u/programs/pDynamo3/installation/shellScripts/environment_bash.com",
                "PROBE\t/home/u/programs/pDynamo3\tsystem\tpython3\tImportError: undefined symbol",
                "PROBE\t/home/u/programs/pDynamo3\tmodule\tpython/3.12.1/gcc.8.5.0\tModuleNotFoundError: No module named 'numpy'",
                "PROBE\t/home/u/programs/pDynamo3\tmodule\tpython/3.10.13/gcc.8.5.0\t3.10.13",
                "EXEC\t/opt/apps/orca/6.1.1/orca\t",
                "NORUN\t/opt/apps/xtb/6.7.0/bin/xtb\terror while loading shared libraries: libmkl.so",
                "EXEC\t/home/u/programs/xtb-6.7.1/bin/xtb\tversion 6.7.1"]) + "\n", ""

    monkeypatch.setattr(disc, "_local_pdynamo_commit", lambda: "abc123")
    r = disc.discover("hpc", transport=Fake())
    p = r["profile"]
    assert p["scheduler"] == "slurm" and p["workdir"] == "/scratch/u/easyhybrid"
    assert p["slurm"]["partition"] == "compute" and p["slurm"]["account"] == "lab"
    assert p["setup"].splitlines() == ["module load python/3.10.13/gcc.8.5.0",
                                       "source ~/programs/pDynamo3/installation/shellScripts/environment_bash.com",
                                       "module load openmpi/4.1.8/gcc.8.5.0/noncuda"]
    assert p["qc_programs"] == {"orca": "/opt/apps/orca/6.1.1/orca", "xtb": "/home/u/programs/xtb-6.7.1/bin/xtb", "dftb": ""}
    assert p["pdynamo_home"] == "~/programs/pDynamo3" and p["sync_pdynamo"] is True
    notes = " ".join(r["notes"])
    assert "outside" not in notes and "bin/ folder" in notes and "needs its module loaded" in notes


def test_ssh_config_hosts(tmp_path):
    from remote.discover import ssh_config_hosts
    cfg = tmp_path / "config"
    cfg.write_text("Host cluster1 hpc2\n  HostName x\nHost *\n  ServerAliveInterval 30\nHost !bad\n")
    assert ssh_config_hosts(str(cfg)) == ["cluster1", "hpc2"]


def test_ssh_config_add_host(tmp_path):
    from remote import sshconfig
    cfg = tmp_path / "config"
    cfg.write_text("Host cluster1\n    HostName cluster1.example.edu\n")
    path, backup = sshconfig.add_host("labserver", "203.0.113.10", "alice", 2222, "~/.ssh/id_ed25519_ueg", path=str(cfg))
    text = cfg.read_text()
    assert text.startswith("Host cluster1\n    HostName cluster1.example.edu\n")
    assert "Host labserver\n    HostName 203.0.113.10\n    User alice\n    Port 2222\n    IdentityFile ~/.ssh/id_ed25519_ueg" in text
    assert oct(os.stat(cfg).st_mode & 0o777) == "0o600" and os.path.exists(backup)
    assert sshconfig.existing_aliases(str(cfg)) == ["cluster1", "labserver"]
    with pytest.raises(ValueError, match="already"):
        sshconfig.add_host("labserver", "h", "u", 22, path=str(cfg))
    with pytest.raises(ValueError, match="port"):
        sshconfig.add_host("x", "h", "u", 70000, path=str(cfg))
    assert "Port" not in sshconfig.block("a", "h", "u", 22)
    cmds = sshconfig.setup_commands("labserver", "203.0.113.10", "alice", 2222, str(tmp_path / "nokey"))
    assert cmds[0].startswith("ssh-keygen") and "-p 2222 alice@203.0.113.10" in cmds[1] and cmds[2] == "ssh labserver hostname"


def test_slurm_memory_check_and_oom_diagnosis(tmp_path):
    """ bigblue, 2026-10-30: --mem=4G for 50 xTB processes -> killed for memory after 13 h,
        reported only as a TypeError of the minimizer. Now refused before submission, and
        a fetched scheduler.log with oom-kill events is explained in the error report. """
    from remote.jobs import RemoteJobs, RemoteError, _slurm_megabytes
    assert _slurm_megabytes("4G") == 4096 and _slurm_megabytes("4000") == 4000 and _slurm_megabytes("") is None
    profile = {"name": "bigblue", "scheduler": "slurm", "slurm": {"mem": "4G"}}
    with pytest.raises(RemoteError) as error:
        RemoteJobs._check_memory(profile, {"processes": 50, "cpus": 50})
    assert "82 MB each" in str(error.value) and "at least 25G" in str(error.value)
    RemoteJobs._check_memory(profile, {"processes": 4, "cpus": 4})                    # 1 GB each: fine
    RemoteJobs._check_memory({"scheduler": "slurm", "slurm": {"mem": ""}}, {"processes": 50})   # partition default
    RemoteJobs._check_memory({"scheduler": "direct"}, {"processes": 50})
    (tmp_path / "scheduler.log").write_text("slurmstepd: error: Detected 9693 oom-kill event(s) in StepId=3521560.batch.\n")
    assert "lack of memory" in RemoteJobs.diagnose(str(tmp_path))
    (tmp_path / "scheduler.log").write_text("all good\n")
    assert RemoteJobs.diagnose(str(tmp_path)) is None
