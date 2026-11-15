# Running simulations on remote machines (ssh, SLURM)

EasyHybrid 3.2 can send its simulations to another machine reached by ssh and bring the results
back. The machine can be a workstation (the job runs as a background process) or a cluster (the
job goes through SLURM). The window stays free while the job runs. When the job ends, its results
land in the output folder you chose, in the same layout a local run would write. They are then
imported as usual, from the *Finished!* dialog or the Process Manager.

The code is in `src/remote/` and has no GTK dependency:

| File | Role |
|---|---|
| `profiles.py` | the machines, saved in `~/.easyhybrid/remotes.json` |
| `transport.py` | ssh and rsync, with one shared connection per host |
| `schedulers.py` | `direct` (nohup) and `slurm` (sbatch/squeue/sacct/scancel) |
| `runner.py` | runs a job on the remote machine |
| `jobs.py` | bundle, submission, status, fetch, cancel, connection test; the registry `~/.easyhybrid/remote_jobs.json` |

## What the remote machine needs

- **ssh access with a key.** EasyHybrid never asks for or stores a password. Use an alias of
  `~/.ssh/config`, which also sets the user, key, port and ProxyJump. Connect once from a
  terminal (`ssh <alias>`) to accept the host key. If the machine is behind a VPN, the VPN must
  be on.
- **pDynamo**, importable by a Python that has numpy, after the setup lines of the profile.
  EasyHybrid copies its own code (`src/pdynamo`, `src/util`, `src/remote`, about 17 MB the first
  time) to `<remote folder>/easyhybrid/src` at each submission; `rsync` sends only what changed.
- **rsync** on both machines.

## Adding a machine: Detect Settings

1. **Simulate → Remote Machines... → New.**
2. Choose the **SSH host**: the aliases of `~/.ssh/config` are listed. For a machine known only by
   address, user and port, use **Add SSH Host...**. It writes an alias to `~/.ssh/config` (with a
   backup; existing aliases are never changed) and shows the commands to type in a terminal:
   - `ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519_<alias>`, to create a key for that machine;
   - `ssh-copy-id -i ~/.ssh/id_ed25519_<alias>.pub -p <port> <user>@<address>`, the only time the
     account's password is asked;
   - `ssh <alias> hostname`, which must answer without a password.

   Turn the VPN on if the machine needs one.
3. Press **Detect Settings**. EasyHybrid asks the machine, in two ssh calls (a few seconds), for:
   - **scheduler:** SLURM if `sbatch` exists, otherwise direct. PBS and LSF are reported as not
     supported yet;
   - **work folder:** the first writable scratch (`$SCRATCH`, `/scratch/$USER`, `/work/$USER`...), with
     `/easyhybrid` added; the home folder only when there is none;
   - **partitions and accounts:** the cluster's default partition and your default account (`sinfo`,
     `sacctmgr`);
   - **Python + pDynamo:** the pDynamo installations in the home folder. Each is tried with the system
     Python, the newest `python/3.x` modules and conda Pythons, and the one that **really imports**
     pDynamo (with numpy and yaml) goes into the setup lines, with the pDynamo environment script.
     The others are reported with the reason they failed;
   - **pDynamo sync:** *Use this computer's pDynamo files* is ticked when the pDynamo there is at the
     same git commit as here;
   - **ORCA, xTB, DFTB+:** from the modules (also when a module's PATH points to the wrong folder),
     the home folder and the PATH. Each candidate is **run** (`--version`), the ones that fail without
     their module's libraries are left out, and the newest version is chosen;
   - **MPI:** with ORCA, a script of yours that sets up ORCA with Open MPI 4.1
     (`~/opt/*env*.sh`, `~/*env*.sh`) is sourced in the setup lines. Otherwise the newest Open MPI 4.1
     module (gcc, non-CUDA) is loaded, for parallel ORCA;
   - **SLURM health:** SLURM is chosen only when its controller answers (`sinfo`). A cluster whose
     controller is down gets the direct scheduler, with a note;
   - **disk:** a note when the work folder has less than 50 GB free.
4. The form is filled and every field can still be changed. The report lists what was found, the
   alternatives (other scratch folders, Pythons, accounts, program versions, MPI modules) and the
   notes. Review it, then **Save** and **Test Connection**.

Nothing is installed or changed on the machine. On the SLURM cluster used for the tests, Detect Settings reproduced the profile that had been set
up by hand (below), including a broken ORCA module and the Python 3.10 module.

The automation layer has the same: `remote_detect host=<alias>` returns the proposal (`save=true`
stores it).

## Setting up a machine: Simulate → Remote Machines...

| Field | Example (a SLURM cluster) | Meaning |
|---|---|---|
| Name | `hpc1` | shown in Simulate → Run Simulations On |
| SSH host | `hpc1` | alias of `~/.ssh/config`, or `user@host` |
| Remote folder | `/scratch/<user>/easyhybrid` | code and jobs; use the scratch file system |
| Python | `python3` | the interpreter, after the setup lines |
| Setup lines | `module load python/3.10.13/gcc.8.5.0`<br>`source ~/programs/pDynamo3/installation/shellScripts/environment_bash.com` | environment of every job; QC programs too (`module load orca/6.1.1`, `PDYNAMO3_ORCACOMMAND=...`) |
| Scheduler | SLURM | or Direct for a workstation |
| CPUs per job | 1 | default; 2D scans and umbrella sampling use them in parallel |
| SLURM | partition `compute`, account `mylab`, time `01:00:00`, memory `4G` | `#SBATCH` options; *Extra* for anything else. The partition list comes from the cluster (`sinfo`, ↻ refreshes it): the default one first, and below it the time limit, the nodes, the idle CPUs and the GPUs of the chosen one |
| Do not fetch | `*.dcd data_collection_traj` | large files left on the remote machine |

**Test Connection** checks these, in order:

1. ssh;
2. the work folder;
3. the SLURM commands;
4. that the code was copied;
5. Python + pDynamo + the EasyHybrid code;
6. a real energy: a small system serialized here is read and computed there, and the result must
   match the local energy;
7. for ORCA and xTB, when the profile gives them, the energy of a water molecule computed by pDynamo
   there with that program (HF/STO-3G, GFN2). It is compared with this computer's when the program
   exists here. This catches a pDynamo that cannot read the program's output, missing libraries and
   wrong paths.

On the test cluster, the system Python (3.6) could not import pDynamo; the `python/3.10.13` module could.

## QC/MM with ORCA, xTB or DFTB+

A system whose QC region uses an external program carries, in its `.pkl`, the scratch folders
of this computer. It does not carry the program itself: pDynamo reads `PDYNAMO3_ORCACOMMAND`,
`PDYNAMO3_XTBCOMMAND` and `PDYNAMO3_DFTBCOMMAND` from the environment when it runs. For a remote
job:

| Issue | What EasyHybrid does |
|---|---|
| Scratch paths of this computer, in the model and its state | The runner moves them to a folder of the job (see [Temporary files of the QC programs](#temporary-files-of-the-qc-programs)). |
| Where the program is | **pDynamo and QC programs on that machine**: the full path of ORCA, xTB and DFTB+ there. The job script exports the `PDYNAMO3_*COMMAND` variable. Leave a field empty if the setup lines already set the variable. |
| Program missing | The job stops at once with "xTB was not found on <node>: give its path ...", not a Python error later. |
| DFTB+ parameters (`skfPath`, a local folder) | The Slater-Koster files are sent with the job and used there. |
| CPUs | The request to SLURM is the parallel processes of the simulation × the threads of the program (xTB `parallel`, ORCA `PALn`). |
| Parallel ORCA | ORCA runs `mpirun`, and Open MPI counts SLURM *tasks*. A `PALn` job therefore asks `--ntasks=n --cpus-per-task=1`. With several processes in the same job, MPI core binding is turned off. The setup lines must load the MPI module ORCA was built with; without `mpirun` the job stops with that advice. |
| A program that fails | The end of the program's own log (pDynamo renames it to `error_*`) is added to the error of the job. |
| pDynamo changed on this computer | **Use this computer's pDynamo Python files and parameters there**: jobs run on a copy of the remote pDynamo (in `<remote folder>/pdynamo3`, compiled modules kept) with the local `.py` files and `parameters/` on top. The installation there is not changed. |

The last option matters. The pDynamo on the test cluster had the same commit as this computer's, but lacked
some local fixes. One of them is the change to `QCModelORCA.py` that reads the orbital energies of
ORCA 6. Without it, ORCA 6.1.1 ran correctly there, but pDynamo got no energy back and the job
failed.

**Modules with a wrong PATH.** On the test cluster, the `orca/6.1.1` module added
`<prefix>/orca/6.1.1/bin` to the PATH, but the binary was `<prefix>/orca/6.1.1/orca`. Give the
real path in the profile. Parallel ORCA there needed the Open MPI 4.1 module (gcc, non-CUDA).

### Temporary files of the QC programs

xTB, ORCA and DFTB+ write their input, output and gradient files in a scratch folder. On the
remote machine this is a folder of the job, never the scratch folder of this computer:

```
<remote folder>/jobs/<job id>/qc_scratch/              one serial job
<remote folder>/jobs/<job id>/qc_scratch/process_<i>/  parallel processes of a 1D scan or US
<remote folder>/jobs/<job id>/qc_scratch/process_<i>_<j>/   2D scans and 2D US
```

For example, `/scratch/<user>/easyhybrid/jobs/20261029-221503-k3x9/qc_scratch/` on a cluster
whose remote folder is `/scratch/<user>/easyhybrid`. The job id is the date and time of the
submission and four random characters, so:

- **Several jobs at once on the same machine, in SLURM or directly, never share files.** Each job has its
  own folder, and each parallel process inside a job has its own subfolder.
- Each job is a process of its own, so the working folder that a program changes into
  affects no other job.
- The scratch files are not fetched: fetching copies only `output/` of the job. They stay on the machine,
  in the job folder, until that folder is cleaned (which deletes the whole job folder there).

**Node-local disk.** To put the files on a faster disk of the compute node, set
`EASYHYBRID_QC_SCRATCH` in the setup lines of the profile, for example
`export EASYHYBRID_QC_SCRATCH=$TMPDIR` (SLURM) or `export EASYHYBRID_QC_SCRATCH=/tmp/$USER`.
The runner then uses `$EASYHYBRID_QC_SCRATCH/<job id>/qc_scratch`, so a fixed path is safe too:
jobs on the same node keep separate folders. EasyHybrid does not delete this folder, and
cleaning the job does not reach it (SLURM usually empties `$TMPDIR` at the end of the job). When a
program fails, its error log is still read there and added to the error of the job.

## Running a simulation remotely

1. **Simulate → Run Simulations On → hpc1 (remote)**. The status bar confirms the choice.
2. Open any simulation window (geometry optimization, scans, MD, NEB, umbrella sampling, normal
   modes, energy refinement, CPR) and press Run as usual. A single-point energy always runs locally.
3. The Process Manager shows the job as `hpc1: sending`, then `queued (slurm 3514228)`, then
   `running`. EasyHybrid asks the machine every 30 s; right-click → *Check Remote Jobs Now* asks
   at once.
4. When the job ends, its results are copied into the output folder. The final structure becomes
   an object, as after a local job, and *Finished!* offers **Import...** for the trajectory.
5. Right-click → *Cancel Remote Job* cancels it (`scancel`, or kill for direct jobs).

Each job gets its own folder `<remote folder>/jobs/<job id>/` with:

- `system.pkl`, `parameters.pkl` and `job.json`;
- the scheduler script;
- `output/`;
- `status.json`, `runner.log` and `scheduler.log`.

The logs are also copied to `<output folder>/.easyhybrid_remote/<job id>/`. For umbrella sampling
from a trajectory or from a multi-frame object, the frames are sent with the job.

### When a job fails

A job that fails (not submitted, an error of the simulation, or ended by SLURM: memory, time
limit, a setup line) opens a **Remote job failed** dialog with the reason; *Show details* has the
full report. The report is also written to `<output folder>/.easyhybrid_remote/<job id>/error.log`
and becomes the log of the job, so **View Log** in the Process Manager opens it. It has:

- the error message (from `sbatch`, the runner or the scheduler state);
- what was asked: CPUs = parallel processes × threads of the QC program, the remote folder;
- the traceback of the runner, when the simulation itself failed;
- the last 60 lines of `runner.log` and `scheduler.log` (for a job ended by SLURM these are
  copied from the machine first; `scheduler.log` usually says why: out of memory, time limit,
  `module: command not found`, ...).

**CPUs per node.** All processes of a job run on one node. Before sending anything, EasyHybrid
compares the CPUs of the job with the largest node of the partition (`sinfo`; `nproc` for direct
jobs). A request that does not fit stops at once with, for example:

```
The job asks for 500 CPUs (50 parallel processes x 10 xTB threads), but the nodes of partition
icomputeq have at most 40. All processes of a job run on one node: lower the number of processes
(Number of CPUs) or the threads of the QC program so that their product is 40 or less.
```

Without this check `sbatch` refuses the job with "CPU count per node can not be satisfied". For
umbrella sampling and 2D scans, the windows are independent, so a common choice is xTB threads = 1
and processes = the CPUs of one node.

## From scripts and AI assistants

The automation layer (`docs/AI_CONTROL.md`) has these commands:

- `remotes`, `remote_save`, `remote_test`;
- `run_on remote=hpc1`;
- `remote=` and `remote_cpus=` on optimize, scan, md, neb and umbrella;
- `remote_jobs`, `remote_fetch`, `remote_cancel`.

`wait` works for remote jobs too.

```bash
ehctl remote_test name=hpc1
ehctl scan --json '{"obj": 1, "remote": "hpc1", "rc1": {"atoms": [5, 17], "step": 0.2, "nsteps": 20}}'
ehctl wait sys=0 && ehctl import_result sys=0 && ehctl pes sys=0
ehctl umbrella sys=0 job=1 force_constant=100 remote=hpc1 remote_cpus=12
```

Without the interface, `remote.jobs.RemoteJobs().submit(system, parameters, "hpc1", folder)`
sends a pDynamo system with a parameter dictionary, and `refresh()` / `fetch()` follow it.

## Tested

On a university SLURM cluster, on 2026-10-28:

- **1D scan of bALA:** 11 s in the queue and running. Barrier 8.29 kJ/mol at 3.66 Å, C5 minimum at
  5.59 kJ/mol, as in the local tutorial run.
- **Umbrella sampling from that scan, 8 CPUs:** PENDING, then COMPLETED. `us_check` OK, minimum
  overlap 0.35.

QC/MM on the same cluster, also on 2026-10-28: bALA with the acetyl group in the QC region and a link
atom at C4–N6. Each job was a geometry optimization; the table compares its first energy there
with the energy here.

| QC program | Remote − local (kJ/mol) | Notes |
|---|---|---|
| ORCA 6.1.1, HF/6-31G* | 2 × 10⁻⁵ | 6.1.0 here |
| ORCA 6.1.1, PAL2 | 2 × 10⁻⁵ | OpenMPI 4.1.8 module |
| xTB 6.7.1, GFN2, 2 threads | 1 × 10⁻⁶ | |
| DFTB+ 23.1, 3ob | 5 × 10⁻⁹ | Slater-Koster files sent with the job |

## Limits (first version)

- **Schedulers:** direct and SLURM only. PBS/Torque, SGE and LSF are new classes in
  `schedulers.py`.
- **pDynamo:** it must already be installed on the machine; an installer is planned. Changes
  made to pDynamo here travel with *Use this computer's pDynamo* only when they are Python
  files or parameters. Changes to compiled modules would need a rebuild there.
- **Jobs from an earlier session:** jobs sent in a session that was closed are followed with
  `remote_jobs` / `remote_fetch` (registry), not yet in the Process Manager.
- **Umbrella windows:** a run uses one node. Array jobs, with one SLURM task per window, are the
  next step for runs larger than a node.
