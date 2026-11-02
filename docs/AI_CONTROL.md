# Controlling EasyHybrid from scripts and AI assistants

EasyHybrid 3.2 can be driven by another program while its window stays open: a Python
script, a shell command, or an AI assistant (Claude Code, Claude Desktop or any client of
the Model Context Protocol). The program sends commands such as "create this system",
"optimize", "scan this coordinate", "give me the energies" or "take a screenshot". It gets
back JSON, and the window shows everything as if the user had done it.

```
 AI assistant ──MCP (stdio)──> mcp_server.py ─┐
 shell ────────────────────────> ehctl.py ─────┼── Unix socket ──> control server ──> AutomationAPI
 Python ───────────────────────> client.py ────┘   (one JSON line      (GTK main thread)     │
                                                    per message)                             ▼
                                                                       p_session / vm_session (same as the menus)
```

All the code is in `src/automation/`:

| File | Role |
|---|---|
| `api.py` | `AutomationAPI`: the commands. Each public method is one command; its signature and docstring are its documentation and its argument schema |
| `server.py` | control server inside EasyHybrid: Unix socket, runs the commands in the GTK main thread, handles modal dialogs |
| `events.py` | event log; in automation mode, messages that dialogs would show are recorded instead |
| `client.py` | client (standard library only), `launch()` |
| `ehctl.py` | command-line client |
| `mcp_server.py` | MCP server (stdio, standard library only) |

## Starting the control server

Use any one of these:

```bash
python3 easyhybrid.py --control                 # dialogs off (see "Dialogs and events")
python3 easyhybrid.py --control-dialogs         # control server, dialogs still shown
EASYHYBRID_CONTROL=1 python3 easyhybrid.py
python3 src/automation/ehctl.py launch          # starts EasyHybrid --control and waits for it
```

In a window that is already open, tick **Edit → AI Control Server**, or type `control on` /
`control off` / `control` (status) in the EasyHybrid terminal. The menu item always shows the
real state, whichever way the server was started.

The socket is `~/.easyhybrid/control.sock`, or `EASYHYBRID_SOCKET` when that variable is set.
A second EasyHybrid running at the same time uses `control-<pid>.sock`. A path longer than the
Unix limit (about 107 bytes) is moved to `/tmp/easyhybrid-<uid>/`. In every case the path in use
is written to `~/.easyhybrid/control.json`, where the clients look for it.

**Security.** The socket has mode 0600 inside a 0700 folder, so only your own user can connect.
The `python` command runs any Python code inside EasyHybrid, with the same rights as you. Do not
expose the socket to other users or forward it over the network.

## Token use and model (Edit → AI Control Server)

Turning the server on from the menu first opens its options. **Edit → AI Control Settings...**
changes them later, also while the server runs. They are saved in `~/.easyhybrid/ai_settings.json`
and also apply to `--control` and to `control on`.

| Option | Economical | Detailed | Effect |
|---|---|---|---|
| Summaries instead of full data tables | on | off | `pes`, `us_check` and `wham` leave out arrays and window tables; `detail="full"` brings them back |
| Simulations return their result in the same call | on | off | optimize/scan/md/neb/umbrella wait (up to *Wait at most*), import and summarize; `wait=false` starts without waiting |
| Short tool descriptions | on | off | the MCP tool list drops from ~7,900 to ~4,500 tokens; the `help` tool gives the full text |
| Compact JSON replies | on | off | no indentation; empty fields left out |
| Screenshot width | 1000 px | full size | a 1920-px window costs ~2,600 tokens, ~700 at 1000 px |

Measured on the runs of the case studies, economical replies are many times smaller:

| Reply | Detailed | Economical |
|---|---|---|
| `pes` of the β-PGM 2D scan | ~3,500 tokens | ~250 tokens |
| `us_check` of 64 windows | ~5,450 tokens | ~150 tokens |
| `wham` | ~830 tokens | ~220 tokens |

**Model.** EasyHybrid cannot change the model of a Claude Code session; it is chosen in Claude
Code itself with `/model`. What the dialog can do is choose the model of the subagent that
operates EasyHybrid:

- With Haiku, Sonnet or Opus chosen, EasyHybrid writes `~/.claude/agents/easyhybrid-operator.md`.
  This is a Claude Code subagent with that model and short operating instructions. The MCP
  instructions then ask the main session to delegate the execution steps to it, and to keep
  planning and interpretation for itself.
- "Same as the Claude Code session" removes the subagent.
- EasyHybrid only overwrites or removes a subagent file it wrote itself.

Both the subagent and the tool descriptions apply to Claude Code sessions started afterwards.

## Clients

### Command line

```bash
alias ehctl="python3 /path/to/EasyHybrid/src/automation/ehctl.py"
ehctl help                      # list of commands
ehctl help scan                 # arguments of one command
ehctl state
ehctl new_system coordinates=examples/tutorials/01_bala/bALA_not_opt.mol force_field=opls name=bALA
ehctl optimize obj=0 rms_gradient=0.05
ehctl wait sys=0
ehctl scan --json '{"obj": 1, "rc1": {"type": "distance", "atoms": [5, 17], "step": 0.2, "nsteps": 15}}'
ehctl terminal command="show sticks obj=1"
ehctl screenshot path=/tmp/view.png
```

Values written as `key=value` are read as JSON when possible (numbers, `true`/`false`, `null`,
`[lists]`, `{dicts}`) and as text otherwise. The output is the JSON reply.

### Python

```python
import sys; sys.path.insert(0, "/path/to/EasyHybrid/src")
from automation.client import EasyHybridClient, launch

launch()                                       # no-op if already running
eh = EasyHybridClient().connect()
system = eh.call("new_system", coordinates="bALA_not_opt.mol", force_field="opls", name="bALA")
job = eh.call("optimize", obj=system["objects"][0])
eh.wait(sys=job["system"], job=job["job"])
opt = eh.call("import_result", sys=job["system"], job=job["job"])
print(eh.call("energy", obj=opt["object"])["energy"])
```

`call()` returns the result or raises `EasyHybridError`. `request()` returns the whole reply,
including the events.

### AI assistants (MCP)

Claude Code:

```bash
claude mcp add easyhybrid -- python3 /path/to/EasyHybrid/src/automation/mcp_server.py
```

Claude Desktop and other MCP clients (`claude_desktop_config.json` or equivalent):

```json
{ "mcpServers": { "easyhybrid": { "command": "python3",
    "args": ["/path/to/EasyHybrid/src/automation/mcp_server.py"] } } }
```

Each command becomes a tool, plus three extra tools:

- `launch` starts EasyHybrid if it is not running.
- `wait` waits for a job to finish.
- `result` collects the result of a pending request.

`screenshot` also returns the image itself, so the assistant can look at the structure or at a
plot. The MCP server reads the tool list from `api.py` without importing GTK, so it starts even
when EasyHybrid is closed.

### Local models (Ollama, LM Studio, llama.cpp, vLLM)

A step-by-step tutorial that combines Claude Code and a local model on a real study:
`docs/site/ai_assistants_tutorial.md` (site: *EasyHybrid with Claude Code and a local model*).

**Edit → AI Assistant (Local Model)...** opens a chat with a model served on this computer or your
network. It operates the window with the same commands as the MCP server, without Claude and
without sending data anywhere.

- **Server:** *Ollama* (`http://localhost:11434`, models listed with `ollama list`) or an
  *OpenAI-compatible* server such as LM Studio, llama.cpp or vLLM (`/v1/chat/completions`). The
  model must support tool calling; Ollama marks such models with "tools" in the list.
- **Reading is free.** `state`, `atoms`, `measure`, `energy`, `pes`, `job_summary`, `us_check`,
  `jobs` and the display commands run at once.
- **Actions are asked first.** Loading files, defining the QC region or fixed atoms, starting
  optimize/scan/md/neb/umbrella, `wham`, `import_result`, `terminal` commands and remote jobs each
  open a Yes/No dialog, unless *Run actions without asking* is ticked.
- **Never offered to a local model:** `python`, `simulate`, the setup of remote machines and
  screenshots (a text model cannot see them).
- **Summaries only.** The tool list is short (about 4,400 tokens, 39 tools), and replies are always
  summaries: local models have small contexts. *Context (tokens)* sets `num_ctx` for Ollama; 16k is
  enough for a few steps.
- The choices are saved in `~/.easyhybrid/ai_settings.json` (`local_*`).

With `gpt-oss:20b` on an RTX 3060 (12 GB, so part of the model on the CPU), the question "barrier,
reaction energy and HOMO–LUMO gap of this 2D xTB scan" took 36 s. The model chose `job_summary` by
itself and gave the right numbers (99.1 kJ/mol, −26.9 kJ/mol, 5.55–7.83 eV). The first question
also loads the model, which can take a minute or more.

Local models of this size answer questions and run short tasks well. Long protocols with
scientific decisions (window spacing and force constant, why a NEB does not converge) are still
better planned with a large model.

The agent is `src/automation/local_agent.py` (standard library only). `LocalAgent(provider, url,
model, execute, confirm).ask(text)` also works from a script.

## Conventions

- **Ids.** `sys` is the id of a pDynamo system and `obj` the id of an object (a set of coordinates
  or a trajectory). Both come from `state`. Without `sys`, a command uses the active system;
  without `obj`, the last object of that system.
- **Atoms.** Atoms are 0-based pDynamo indices, the same numbers as in the reaction coordinate
  boxes. Give them as `[0, 4, 7]` or `"0-10,15"`. `atoms` lists them with name, residue, chain,
  element, coordinates, charge and QC/fixed flags. It accepts the filters `chain`, `resi`, `resn`,
  `name` and `symbol`, plus `within` + `of` for a distance cut-off.
- **Frames.** `frame` is 0-based; when omitted, the last frame is used.
- **Units.** Energies are in kJ/mol, distances in Å, angles in degrees.
- **Jobs.** Simulations return a job (`system`, `job`) as soon as they start:
  1. call `wait` until the job finishes;
  2. call `import_result` to load the trajectory or scan as a new object;
  3. call `pes` to get the energies of a scan, NEB or energy refinement as numbers.

  A system runs one job at a time, as in the interface.
- **Output folder.** Without `folder`, results go to the system's working folder.

## Dialogs and events

In automation mode (`--control`, the default) these messages are not shown as dialogs:

- information and error dialogs;
- "QC Model defined successfully";
- the "Finished!" dialog of a job.

Each message is recorded as an **event** instead, and the reply to the command that caused it
includes it under `"events"`. `events since=N` lists all events, for example
`kinds=["finished","job_error"]`.

Question dialogs (Yes/No) and other windows that wait for an answer are still shown. When a
command opens one, the reply comes back after about a second as
`{"ok": true, "pending": "r3", "dialogs": [...]}`. Each dialog is listed with its text and
buttons. To continue:

1. answer it with `dialog_respond button="Yes"` (or `response=`, or `text=` for a text field);
2. get the result of the original command with `result request=r3`.

## Commands

| Command | What it does | Arguments |
|---|---|---|
| `state` | Systems, objects, running jobs and open dialogs of the session. | |
| `load` | Opens a file as File > Open (.easy session, .pkl/.yaml system, PDB/XYZ/MOL2...). | path |
| `new_system` | Creates a system as File > New (force field and coordinates). | coordinates, force_field, name, tag, topology, parameters, parameter_folder |
| `save_session` | Saves the session (.easy). | path |
| `activate` | Makes a system the active one. | sys |
| `atoms` | Lists atoms (index, name, residue, chain, element, xyz, charge, QC/fixed). | obj, sys, atoms, chain, resi, resn, name, symbol, within, of, by_residue, frame, limit |
| `measure` | Distance (2 atoms), angle (3) or dihedral (4); one frame or every frame. | atoms, obj, sys, frame, frames |
| `set_qc` | Defines the QC region and the QC method (pDynamo semiempirical, ORCA, xTB). | atoms, method, charge, multiplicity, engine, filters, orca_keywords, gfn, cpus |
| `fix_atoms` | Sets the fixed atoms, e.g. everything farther than 10 Å from the active site. | atoms, free_within, of, by_residue, clear |
| `energy` | Single-point energy with its terms. | sys, obj, frame |
| `optimize` | Geometry optimization. | obj, frame, optimizer, max_iterations, rms_gradient, name, folder, save_trajectory |
| `scan` | Relaxed surface scan, 1D or 2D. | rc1, rc2, obj, frame, name, folder, cpus, optimizer, max_iterations, rms_gradient |
| `md` | Molecular dynamics. | obj, steps, time_step, integrator, temperature, collision_frequency, name, save_frequency, ... |
| `neb` | Nudged elastic band between two objects. | reactant, product, images, spring, spline, name, cpus, ... |
| `umbrella` | Umbrella sampling: windows from a scan trajectory (parallel) or from one structure (serial). | source or sys+job, obj, rc1, rc2, force_constant, cpus, steps_eq, steps_dc, temperature, optimize, name |
| `us_check` | Diagnostics of a run (as the Umbrella Sampling Browser computes them): windows that slid from their targets, overlap gaps, convergence by halves, and a verdict. | folder or sys+job, halves |
| `wham` | PMF with pDynamo's WHAM, as numbers: barrier (from the first window), minima, reaction free energy; 2D: the lowest-barrier path. | folder or sys+job, bins, bins_rc2, temperature, resamples, path_from, path_to |
| `us_browser` | Opens the Umbrella Sampling Browser on a run (then `screenshot target="Umbrella Sampling Browser"`). | folder or sys+job |
| `dock` | Docking with dyffdock (pDynamo + DYFF), as Extras → Run DYFFDock: ligands from mol2 files or atoms of the system, rigid or ensemble receptor; waits and returns the best poses (scores in kJ/mol, contacts). | ligands, ligand_atoms, sys, obj, frame, frames, center, size, pocket, pocket_margin, search, runs, seed, options, folder, show_poses, wait |
| `pockets` | Pockets of a receptor system, largest first (the cavities of Analysis → Cavities and Pockets, of the receptor docked into), each with its lining residues, its overlap with given ligands and its docking box; `dock pocket=n` uses that box. | sys, obj, frame, leave_out, ligand_atoms, ligands, pocket_margin, detail |
| `dock_results` | Results of a dyffdock run folder; `show_poses` loads each ligand with its poses as frames. | folder, show_poses, detail |
| `simulate` | Any simulation from its full parameter dictionary (normal modes, energy refinement, ...). | parameters, obj, frame |
| `remote_detect` | Asks a machine (ssh) for scheduler, partitions, accounts, scratch, Python + pDynamo, QC programs; proposes a profile. | host, name, save |
| `remotes` / `remote_save` / `remote_test` | Remote machines (ssh, SLURM): list, add, test (see docs/REMOTE_JOBS.md). | name, host, workdir, setup, scheduler, partition, account, ... |
| `run_on` | Where the next simulations run: this computer or a remote machine. | remote |
| `remote_jobs` / `remote_fetch` / `remote_cancel` | Remote jobs: state, results now, cancel. Simulations also take `remote=` and `remote_cpus=`. | refresh, remote_job |
| `jobs` / `job_status` | Jobs of a system / one job (status, log file, error, resulting object). | sys, job |
| `wait` | Waits for a job to finish (server side). | sys, job, timeout |
| `import_result` | Imports the trajectory of a finished job as a new object. | sys, job, name |
| `job_summary` | What the log of any job says, as numbers (the "Plots" tab of View Log): convergence, MD averages, barriers, NEB top, umbrella verdict, imaginary frequencies, ORCA/xTB energies and gaps. | sys, job, logfile, detail |
| `pes` | Energy profile or surface of a scan/NEB/refinement as numbers: barriers (local maxima), minima; 2D: lowest-barrier path. | sys, job, logfile, relative, path_from, path_to |
| `terminal` | Any EasyHybrid terminal command: show, hide, color, center, select, zoom, dyff, opls, solvate, namd... | command |
| `set_visible` | Shows or hides an object (`only=true` hides the others). | obj, visible, only |
| `set_frame` | Shows a frame. | frame |
| `screenshot` | PNG of the 3D view, of the main window or of any window by title. | path, target, scale |
| `windows` | Titles of the open windows. | |
| `dialogs` / `dialog_respond` | Open dialogs; answer one. | button, dialog, response, text |
| `events` | Recorded events. | since, kinds |
| `python` | Runs Python inside EasyHybrid (`main`, `p_session`, `vm_session`, `api`, `np`; assign `result`). | code |
| `ping` | Connection check, version. | |

A reaction coordinate for `scan` is a dict:

```json
{"type": "distance" | "multiple_distance" | "multiple_distance_4" | "dihedral",
 "atoms": [i, j], "step": 0.1, "nsteps": 20, "force_constant": 4000,
 "mass_weighted": false, "start": null}
```

- `multiple_distance` (3 atoms) is d(i–j) − d(j–k): atom j moves from i to k.
- `multiple_distance_4` (4 atoms) is d(i–j) − d(k–l).
- `start` defaults to the value in the starting coordinates.
- With `rc2` the scan is 2D and runs on `cpus` processes.

## Example: the bALA scan of the tutorial, by commands

```bash
ehctl new_system coordinates=examples/tutorials/01_bala/bALA_not_opt.mol force_field=opls name=bALA tag=bALA
ehctl optimize obj=0 folder=$PWD/work
ehctl wait sys=0
ehctl import_result sys=0                      # -> object 1 (single structure)
ehctl atoms obj=1 name=O5,H17                  # the C7eq hydrogen bond: atoms 5 and 17
ehctl scan --json '{"obj": 1, "folder": "'$PWD'/work", "name": "bALA_scan",
                    "rc1": {"atoms": [5, 17], "step": 0.2, "nsteps": 20}}'
ehctl wait sys=0
ehctl import_result sys=0 name=bALA_scan       # -> object 3, 20 frames
ehctl pes sys=0                                # energies, barrier (local_maxima), minima
ehctl measure obj=3 atoms='[5,17]' frames=all
ehctl screenshot target="main"
```

## Example: free energy by umbrella sampling

Umbrella sampling continues from the scan above (job 1 of system 0). The windows are its 20
frames, and the reaction coordinate is read from the scan's log:

```bash
ehctl umbrella sys=0 job=1 force_constant=100 cpus=12 name=bALA_US   # -> job 2
ehctl wait sys=0 job=2
ehctl us_check sys=0 job=2       # verdict OK / warning / fail, with the reasons
ehctl wham sys=0 job=2           # barrier_from_first_window, local_minima, reaction_free_energy
ehctl us_browser sys=0 job=2 && ehctl screenshot target="Umbrella Sampling Browser"
```

`us_check` applies the criteria of the Browser:

- **fail** when a window's mean lies more than 3 σ from its target (the restraint did not hold the
  coordinate) or when neighbouring windows overlap by less than 0.03;
- **warning** when the PMFs of the two halves of the data differ by more than 4 kJ/mol between
  the first and the last target.

Spacing and force constant go together: windows d apart overlap when d is about 2–3 σ, with
σ = √(kT/2k). For example, k = 600 kJ mol⁻¹ Å⁻² gives σ = 0.046 Å at 300 K. On the chorismate
mutase runs of the case study, `us_check` fails the first protocol (windows 16–19 slid, halves
10.5 kJ/mol apart) and passes the corrected one. `wham` gives barriers of 121.8 and 120.8 kJ/mol
(the case study reports about 122 and 120).

2D surfaces (`pes` of a 2D scan, `wham` of a 2D run) use cells `[i, j]` = (RC1 step, RC2 step)
and report the lowest-barrier path between the two minima. Isolated artefact cells can attract
that path: then give `path_from` and `path_to`. On the β-PGM scan,
`path_from=[2,2] path_to=[13,18]` gives the 22.5 kJ/mol barrier of the case study.

## Adding a command

Add a method to `AutomationAPI` in `src/automation/api.py`:

- Use keyword arguments with defaults, and return JSON-compatible values.
- Raise `APIError` for errors the user can fix.
- Write a one-line docstring summary and an `Args:` block.

The method then appears in `ehctl help`, in `describe` and in the MCP tool list, with no further
registration. It runs in the GTK main thread, so it can use the interface; long work belongs in
a job (`run_simulation`) or a thread.

## Limits

- Normal modes, energy refinement, the transition-state search and CPR have no short command
  yet. Run them with `simulate` and the parameter dictionary their window builds (**Export
  Script** shows it).
- Most analysis windows (RMSD, RDF, trajectory analysis...) have no command. Open them with
  `python` and capture them with `screenshot target="<title>"`, or analyze the files directly.
  `pes`, `us_check` and `wham` cover scans, paths and umbrella sampling.
- `us_check` and `wham` read files and run in a worker thread, so the window stays responsive.
  The other commands run in the interface thread.
- The Builder is not exposed.
- The Unix socket needs Linux or macOS.
