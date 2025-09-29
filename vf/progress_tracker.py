# vf/progress_tracker.py
import json, os, time
from typing import Dict, Union, Iterable

# Parent folder where each run will keep its progress.json
BASE_OUT = os.path.join("media", "vf_uploads")

def _run_dir(run_id: str) -> str:
    d = os.path.join(BASE_OUT, run_id)
    os.makedirs(d, exist_ok=True)
    return d

def _progress_path(run_id: str) -> str:
    return os.path.join(_run_dir(run_id), "progress.json")

def _write_atomic(path: str, data: dict) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f)
    os.replace(tmp, path)  # atomic on the same filesystem

def init_progress(run_id: str) -> None:
    """Create/reset progress.json for a new run."""
    data = {"percent": 0, "stage": "Queued", "log": [], "done": False, "error": None, "ts": time.time()}
    _write_atomic(_progress_path(run_id), data)

def read_progress(run_id: str) -> Dict:
    """Read current progress safely (returns defaults if missing)."""
    path = _progress_path(run_id)
    if not os.path.exists(path):
        return {"percent": 0, "stage": "Queued", "log": [], "done": False, "error": None, "ts": time.time()}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def update_progress(
    run_id: str,
    *,
    percent: Union[int, None] = None,
    stage: Union[str, None] = None,
    log: Union[str, Iterable[str], None] = None,
    done: Union[bool, None] = None,
    error: Union[str, None] = None,
) -> None:
    """Update fields in progress.json (append logs instead of overwriting)."""
    data = read_progress(run_id)
    if percent is not None: data["percent"] = int(percent)
    if stage   is not None: data["stage"]   = stage
    if log:
        if isinstance(log, (list, tuple, set)):
            data["log"].extend([str(x) for x in log])
        else:
            data["log"].append(str(log))
    if done  is not None: data["done"]  = bool(done)
    if error is not None: data["error"] = str(error)
    data["ts"] = time.time()
    _write_atomic(_progress_path(run_id), data)
