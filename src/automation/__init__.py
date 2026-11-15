#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Automation package (2026-10-28): lets an external program -- a script, a
#  terminal command or an AI agent -- drive a running EasyHybrid.
#
#    events.py      event log + "automation mode" flag (no GTK)
#    api.py         AutomationAPI: the commands, JSON in / JSON out
#    server.py      control server: Unix socket inside the GUI process
#    client.py      socket client (standard library only)
#    ehctl.py       command-line client:  python3 src/automation/ehctl.py state
#    mcp_server.py  Model Context Protocol server (stdio) for AI assistants
#
#  See docs/AI_CONTROL.md.
