#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Remote jobs (2026-10-28): runs the simulations of EasyHybrid on other
#  machines reached by ssh -- a workstation, or a cluster through its
#  scheduler. GTK-free, so the same code serves the interface, the
#  automation layer and scripts.
#
#    profiles.py    the machines (~/.easyhybrid/remotes.json)
#    transport.py   ssh / rsync with a shared connection (ControlMaster)
#    schedulers.py  how a job is started and followed: direct (nohup), SLURM
#    runner.py      runs a job on the remote machine (pDynamo + src/pdynamo)
#    jobs.py        job bundle, submission, status, fetch, registry
#
#  On the remote machine: <workdir>/easyhybrid/src (a copy of src/pdynamo,
#  src/util and src/remote, updated at each submission) and
#  <workdir>/jobs/<job id>/ (system.pkl, parameters.pkl, job.json, the
#  scheduler script, output/, status.json). See docs/REMOTE_JOBS.md.
