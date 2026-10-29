#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  Remote jobs in the session (2026-10-28): run_simulation() sends the job to
#  another machine (src/remote) when one is chosen in Simulate > Run
#  Simulations On, or with parameters["remote"]. The job gets a row in the
#  Process Manager and an entry in the system's job history like a local
#  one; when it ends, its results are fetched into the output folder and
#  the usual "Finished!" (with Import...) follows.
#
#  Network work (transfer, status, fetch) runs in threads; the interface is
#  only touched from GLib idle callbacks.

import os
import threading

from gi.repository import GLib

from util.debug import dprint

# . simulation types that may run remotely (a single-point energy stays local)
REMOTE_TYPES = ("Geometry_Optimization", "Relaxed_Surface_Scan", "Advanced_Relaxed_Surface_Scan",
                "Molecular_Dynamics", "Nudged_Elastic_Band", "Umbrella_Sampling", "Normal_Modes",
                "Energy_Refinement", "Conjugate_Peak_Refinement")
POLL_SECONDS = 30


class pRemote:
    """ Mixin of pDynamoSession. """

    remote_target = None            # . profile name chosen in Simulate > Run Simulations On

    def _remote_manager(self):
        if getattr(self, "_remote_jobs_manager", None) is None:
            from remote.jobs import RemoteJobs
            self._remote_jobs_manager = RemoteJobs()
            self._remote_rows = {}
        return self._remote_jobs_manager

    # ------------------------------------------------------------------ submit
    def run_remote(self, parameters, profile_name, cpus=None):
        """ Prepares the job here and sends it in a thread. Returns the job
            record (state "preparing"); False when nothing was started. """
        manager = self._remote_manager()
        system = self.psystem.get(self.active_id)
        if system is None:
            self.main.simple_dialog.error(msg="No active system selected. Please load or select a system first.")
            return False
        folder = parameters.get("folder") or getattr(system, "e_working_folder", None) or os.getcwd()
        if not os.path.isdir(folder):
            self.main.simple_dialog.error(msg="The output folder does not exist:\n{}".format(folder))
            return False
        parameters["folder"] = folder
        step = system.e_step_counter
        key6 = self.vm_session.gen_random_tag_string(length=6)
        sim_type = parameters["simulation_type"]
        backup = {k: v for k, v in parameters.items() if k not in ("system", "queue")}
        system.e_job_history[step] = {
            "backup_parameters": backup, "new_vobject": None, "energy": None, "e_id": system.e_id,
            "simulation_type": sim_type, "error": None, "step_counter": step, "started": None, "ended": None,
            "status": "Queued", "potential": self.get_hamiltonian_type(system.e_id), "key6": key6,
            "remote": profile_name, "logfile": None,
        }
        parameters["system"] = system
        parameters = self._apply_restraints(parameters)          # . the restraint model travels in system.pkl
        parameters.pop("system", None)
        try:
            job = manager.prepare(system, parameters, profile_name, folder, cpus=cpus)
        except Exception as error:
            del system.e_job_history[step]
            self.main.simple_dialog.error(msg="Could not prepare the remote job:\n{}".format(error))
            return False
        system.e_job_history[step]["remote_job"] = job["id"]
        self.main.process_manager_window.open_window()
        treeiter = self.main.process_manager_window.add_new_process(
            system=system, _type=sim_type, potential=system.e_job_history[step]["potential"],
            status="{}: sending".format(profile_name))
        self._remote_rows[job["id"]] = (system.e_id, step, treeiter)
        self.main.bottom_notebook.status_teeview_add_new_item(
            message="{} {} - sent to {}".format(sim_type, step, profile_name), system=system)
        # . the next job of this system gets the next number, as after a local job
        system.e_step_counter += 1
        self.changed = True

        def work():
            try:
                record = manager.send(job["id"])
                GLib.idle_add(self._remote_show, record)
            except Exception as error:
                GLib.idle_add(self._remote_failed, job["id"], "Not submitted: {}".format(error))
            self._remote_start_polling()

        threading.Thread(target=work, name="easyhybrid-remote-send", daemon=True).start()
        return job

    # ------------------------------------------------------------------ status
    def _remote_start_polling(self):
        if getattr(self, "_remote_poll_source", None) is None:
            self._remote_poll_source = GLib.timeout_add_seconds(POLL_SECONDS, self._remote_poll_tick)

    def _remote_poll_tick(self):
        ids = [i for i in self._remote_rows if not self._remote_is_final(i)]
        if not ids:
            self._remote_poll_source = None
            return False
        threading.Thread(target=self._remote_poll_work, args=(ids,), daemon=True).start()
        return True

    def _remote_is_final(self, job_id):
        from remote.jobs import FINAL_STATES
        try:
            return self._remote_manager().get(job_id).get("state") in FINAL_STATES and \
                job_id not in getattr(self, "_remote_pending_fetch", set())
        except Exception:
            return True

    def remote_refresh_now(self, job_ids=None):
        """ Status check now (Process Manager / automation), in a thread. """
        ids = job_ids or list(getattr(self, "_remote_rows", {}))
        if ids:
            threading.Thread(target=self._remote_poll_work, args=(ids,), daemon=True).start()

    def _remote_poll_work(self, ids):
        manager = self._remote_manager()
        try:
            records = manager.refresh(ids)
        except Exception as error:
            dprint("remote status:", error)
            return
        for record in records:
            GLib.idle_add(self._remote_show, record)
            if record.get("state") in ("finished", "error") and not record.get("fetched"):
                pending = self.__dict__.setdefault("_remote_pending_fetch", set())
                if record["id"] in pending:
                    continue
                pending.add(record["id"])
                try:
                    fetched = manager.fetch(record["id"])
                    GLib.idle_add(self._remote_done, fetched)
                except Exception as error:
                    GLib.idle_add(self._remote_failed, record["id"], "Could not fetch the results: {}".format(error))
                finally:
                    pending.discard(record["id"])
            elif record.get("state") in ("failed", "cancelled"):
                GLib.idle_add(self._remote_failed, record["id"],
                              record.get("error") or "The job ended without results ({}).".format(
                                  record.get("scheduler_state") or record.get("state")))

    def _remote_row(self, job_id):
        return getattr(self, "_remote_rows", {}).get(job_id)

    def _remote_show(self, record):
        row = self._remote_row(record["id"])
        if row is None:
            return False
        e_id, step, treeiter = row
        state = record.get("state")
        if state in ("queued", "running", "preparing", "submitted"):
            label = "{}: {}".format(record["remote"], state)
            if record.get("scheduler_id"):
                label += " ({} {})".format(record.get("scheduler"), record["scheduler_id"])
            self.main.job_history_liststore[treeiter][5] = label
            system = self.psystem.get(e_id)
            if system is not None and step in system.e_job_history:
                system.e_job_history[step]["status"] = "Running" if state == "running" else "Queued"
                system.e_job_history[step]["remote_state"] = state
        return False

    def _remote_done(self, record):
        """ Fetched: the job history, the final coordinates as an object, Finished!. """
        row = self._remote_row(record["id"])
        if row is None:
            return False
        e_id, step, treeiter = row
        system = self.psystem.get(e_id)
        if system is None:
            return False
        if record.get("state") == "error":
            return self._remote_failed(record["id"], record.get("error") or "The remote job failed.")
        entry = system.e_job_history.setdefault(step, {})
        entry.update({"logfile": record.get("logfile"), "warning": record.get("warning"), "status": "Finished",
                      "remote_state": "finished"})
        final = record.get("final_coordinates")
        if final and os.path.isfile(final):
            try:
                from pCore import Unpickle
                system.coordinates3 = Unpickle(final)
                name = "{}_{}".format(entry.get("simulation_type", "remote"), step)
                vobject = self._add_vismol_object_to_easyhybrid_session(system=system, name=name, key6=entry.get("key6"))
                entry["new_vobject"] = True
                vobject.results = {"remote_job": record["id"]}
            except Exception as error:
                dprint("remote final coordinates:", error)
        pm = self.main.process_manager_window
        pm.set_time(treeiter, False, True)
        pm.set_status(treeiter, "Finished", warning=record.get("warning"))
        self.main.bottom_notebook.status_teeview_add_new_item(
            message="{} {} - finished on {}, results in {}".format(entry.get("simulation_type"), step,
                                                                   record["remote"], record["local_folder"]),
            system=system)
        self.changed = True
        return False

    def _remote_failed(self, job_id, message):
        row = self._remote_row(job_id)
        try:
            self._remote_manager()._update(job_id, state="error" if "Not submitted" in message else
                                           self._remote_manager().get(job_id).get("state"), error=message)
        except Exception:
            pass
        from automation import events as automation_events
        if row is None:
            automation_events.record("job_error", message, title="Remote job", remote_job=job_id)
            return False
        e_id, step, treeiter = row
        system = self.psystem.get(e_id)
        if system is not None and step in system.e_job_history:
            system.e_job_history[step].update({"status": "Error!", "error": message})
        self.main.process_manager_window.set_status(treeiter, "Error!")
        automation_events.record("job_error", message, title="Remote job", system=e_id, job=step, remote_job=job_id)
        self.main.bottom_notebook.status_teeview_add_new_item(message="Remote job {}: {}".format(job_id, message),
                                                              system=system)
        return False

    # ------------------------------------------------------------------ actions
    def remote_cancel(self, job_id):
        manager = self._remote_manager()

        def work():
            try:
                manager.cancel(job_id)
                GLib.idle_add(self._remote_failed, job_id, "Cancelled.")
            except Exception as error:
                GLib.idle_add(self._remote_failed, job_id, "Could not cancel: {}".format(error))

        threading.Thread(target=work, daemon=True).start()

    def remote_job_of_row(self, e_id, step):
        system = self.psystem.get(e_id)
        if system is None:
            return None
        return (system.e_job_history.get(step) or {}).get("remote_job")
