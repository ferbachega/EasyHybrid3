# Notes for AI agents

## Operating EasyHybrid (running simulations, inspecting results)

Do not script the GUI with clicks or screen grabs. EasyHybrid has a control server: start it
with `python3 easyhybrid.py --control` or `python3 src/automation/ehctl.py launch`, then send
commands with the `ehctl` client (or the Python client, or the MCP server). Full guide:
[docs/AI_CONTROL.md](docs/AI_CONTROL.md).

```bash
EH="python3 src/automation/ehctl.py"
$EH state                    # systems (sys=), objects (obj=), jobs, open dialogs
$EH help                     # all commands; "$EH help scan" for the arguments of one
$EH new_system coordinates=examples/tutorials/01_bala/bALA_not_opt.mol force_field=opls name=bALA
$EH optimize obj=0 && $EH wait sys=0 && $EH import_result sys=0
$EH screenshot target=main   # then read the PNG to look at it
```

- Atom indices are 0-based. Energies are in kJ/mol, distances in Å.
- Simulations return a job immediately: call `wait`, then `import_result`, then `pes` for the
  energies as numbers.
- For umbrella sampling, call `umbrella` (from a scan job), then `wait`, then `us_check` (do not
  trust the PMF when the verdict is "fail"), then `wham`.
- Remote machines (ssh, SLURM): `remotes`, `remote_test`, then `remote="<name>"` on a simulation
  command; `wait` works, results come back to the output folder (docs/REMOTE_JOBS.md).
- When a reply says `"pending"`, a dialog is open: call `dialogs`, then `dialog_respond`, then
  `result`.
- Write "QC/MM", never "QM/MM".

## Working on the code

- `src/automation/api.py`: each public method of `AutomationAPI` is one command. Its signature
  and docstring (`Args:` block) generate the CLI help and the MCP schema, so new features should
  get a method there.
- Tests: `python3 -m pytest tests -q --ignore=tests/test_bond_order_matching.py`.
- `src/graphics_engine` is a separate, nested git repository.
- Every string the user sees in the interface is in English.
