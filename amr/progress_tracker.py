# amr/progress_tracker.py
import json, os, time
from typing import Dict, Union, Iterable, Optional

# Base folder: MEDIA_ROOT/amr_runs/<RUN_ID>/progress.json
MEDIA_ROOT = os.environ.get("MEDIA_ROOT", "media")
BASE_OUT   = os.path.join(MEDIA_ROOT, "amr_runs")

def _run_dir(run_id: str) -> str:
    d = os.path.join(BASE_OUT, str(run_id))
    os.makedirs(d, exist_ok=True)
    return d

def progress_path(run_id: str) -> str:
    """Absolute path to the progress.json for a run."""
    return os.path.join(_run_dir(run_id), "progress.json")

def _write_atomic(path: str, data: dict) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f)
    os.replace(tmp, path)

def _default_payload() -> Dict:
    # Shape expected by progress.html
    return {
        "percent": 0,
        "message": "Queued",
        "done": False,
        "error": None,
        "redirect": None,
        "ts": time.time(),
        # optional debugging field (unused by UI but handy)
        "log": [],
    }

def init_progress(run_id: str) -> None:
    """Create/reset progress.json for a new AMR run."""
    _write_atomic(progress_path(run_id), _default_payload())

def read_progress(run_id: str) -> Dict:
    """Read current progress safely (returns defaults if missing/corrupt)."""
    path = progress_path(run_id)
    data = _default_payload()
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                loaded = json.load(f)
                if isinstance(loaded, dict):
                    data.update(loaded)
        except Exception:
            # Keep defaults on any read/parse error
            pass
    return data

def update_progress(
    run_id: str,
    *,
    percent: Optional[int] = None,
    message: Optional[str] = None,   # <- what the UI displays
    done: Optional[bool] = None,
    error: Optional[str] = None,
    redirect: Optional[str] = None,  # optional final URL
    log: Union[str, Iterable[str], None] = None,  # optional debug log
) -> None:
    """Merge-update fields in progress.json (append logs instead of overwriting)."""
    data = read_progress(run_id)
    if percent is not None:  data["percent"]  = int(percent)
    if message is not None:  data["message"]  = str(message)
    if done is not None:     data["done"]     = bool(done)
    if error is not None:    data["error"]    = str(error)
    if redirect is not None: data["redirect"] = str(redirect)
    if log:
        if isinstance(log, (list, tuple, set)):
            data["log"].extend([str(x) for x in log])
        else:
            data["log"].append(str(log))
    data["ts"] = time.time()
    _write_atomic(progress_path(run_id), data)

def env_for(run_id: str) -> Dict[str, str]:
    """
    Environment mapping to launch AMR_pred.py so it can write progress:
    - RUN_ID: used by the script
    - MEDIA_ROOT: base media path for this Django instance
    """
    env = os.environ.copy()
    env["RUN_ID"] = str(run_id)
    env["MEDIA_ROOT"] = str(MEDIA_ROOT)
    return env
