import io
import os
import sys
import uuid
import shutil
import zipfile
import subprocess
import json
import time
import threading
import csv

from django.conf import settings
from django.http import FileResponse, Http404, JsonResponse
from django.shortcuts import render, redirect
from django.urls import reverse

from .forms import UploadFastaForm


# ============================= Helpers =====================================

def _bytes_to_human(n: int) -> str:
    for unit in ("bytes", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "bytes" else f"{n:.1f} {unit}"
        n /= 1024


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def clean_dir(path: str) -> None:
    if not os.path.isdir(path):
        return
    for name in os.listdir(path):
        p = os.path.join(path, name)
        try:
            if os.path.isfile(p):
                os.remove(p)
        except Exception:
            pass


def copy_with_retry(src: str, dst: str, attempts: int = 40, delay: float = 0.5) -> None:
    if os.path.abspath(src) == os.path.abspath(dst):
        return
    last_err = None
    for _ in range(attempts):
        try:
            shutil.copy2(src, dst)
            return
        except (PermissionError, OSError) as e:
            last_err = e
            time.sleep(delay)
    if last_err:
        raise last_err


def _write_progress(run_dir, pct, msg, done=False, error=None, redirect=None,
                    attempts: int = 40, delay: float = 0.15):
    data = {"percent": int(pct), "message": msg, "done": done, "error": error, "redirect": redirect}
    os.makedirs(run_dir, exist_ok=True)
    path = os.path.join(run_dir, "progress.json")
    tmp = path + ".tmp"
    for i in range(attempts):
        try:
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(data, f)
            os.replace(tmp, path)
            return
        except (PermissionError, OSError):
            time.sleep(delay)
            if i % 5 == 4:
                delay = min(delay * 1.4, 1.0)
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f)
    except Exception:
        pass


def _first_existing(*paths: str):
    for p in paths:
        if p and os.path.isfile(p):
            return p
    return None


# ---- FASTA subset from multiclass CSV IDs ---------------------------------

def _ids_from_csv(csv_path: str, id_col: str = "ID") -> set[str]:
    ids = set()
    if not os.path.isfile(csv_path):
        return ids
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = (row.get(id_col) or "").strip()
            if not val:
                continue
            ids.add(val.split()[0])
    return ids


def _write_subset_fasta(src_fasta: str, keep_ids: set[str], out_fasta: str) -> int:
    if not keep_ids or not os.path.isfile(src_fasta):
        return 0
    written = 0
    with open(src_fasta, "r", encoding="utf-8") as fi, \
         open(out_fasta, "w", encoding="utf-8") as fo:

        cur_header, cur_id, cur_seq = None, None, []

        def flush():
            nonlocal written
            if cur_header and cur_id in keep_ids:
                fo.write(cur_header)
                fo.writelines(cur_seq)
                if not cur_seq or not cur_seq[-1].endswith("\n"):
                    fo.write("\n")
                written += 1

        for line in fi:
            if line.startswith(">"):
                if cur_header is not None:
                    flush()
                cur_header = line
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_header is not None:
            flush()

    if written == 0 and os.path.isfile(out_fasta):
        try:
            os.remove(out_fasta)
        except Exception:
            pass

    return written


# ========================== Background job =================================

BIN_NAMES_IN   = {"binary_prediction.csv", "prediction_vf_binary.csv"}
MULTI_NAMES_IN = {"multiclass_only_predictions.csv", "prediction_vf_multiclass.csv"}

BIN_NAME_OUT   = "binary_prediction.csv"
MULTI_NAME_OUT = "multiclass_only_predictions.csv"

def _run_vf_job(inp_path: str, run_id: str):
    run_dir = ensure_dir(os.path.join(settings.MEDIA_ROOT, "vf_runs", run_id))
    public_features = ensure_dir(os.path.join(settings.MEDIA_ROOT, "features_vf"))
    public_pred     = ensure_dir(os.path.join(settings.MEDIA_ROOT, "pred_vf"))

    try:
        _write_progress(run_dir, 10, "Queued…"); time.sleep(0.15)
        _write_progress(run_dir, 25, "Extracting peptide features…")
        _write_progress(run_dir, 40, "Predicting PTMs…")
        _write_progress(run_dir, 55, "Computing ESM2 embeddings…")

        script_path = os.path.join(settings.BASE_DIR, "vf", "VF_pred.py")
        if not os.path.isfile(script_path):
            raise FileNotFoundError(f"Pipeline script not found: {script_path}")

        env = os.environ.copy()
        env["RUN_ID"] = str(run_id)
        env["MEDIA_ROOT"] = str(settings.MEDIA_ROOT)

        subprocess.run([sys.executable, script_path, inp_path],
                       check=True, cwd=settings.BASE_DIR, env=env)

        pref_features_dir   = os.path.join(settings.MEDIA_ROOT, "features_vf")
        legacy_features_dir = os.path.join(settings.MEDIA_ROOT, "vf_uploads", "features_vf")
        src_features_dir = pref_features_dir if (os.path.isdir(pref_features_dir) and os.listdir(pref_features_dir)) else legacy_features_dir

        pred_dir_primary = os.path.join(settings.MEDIA_ROOT, "pred_vf")
        pred_dir_legacy  = os.path.join(settings.MEDIA_ROOT, "vf_uploads", "pred_vf")
        predict_vf_dir   = os.path.join(settings.MEDIA_ROOT, "predict_vf")

        def find_any(dirpath, names):
            if not os.path.isdir(dirpath):
                return None
            for n in names:
                p = os.path.join(dirpath, n)
                if os.path.isfile(p):
                    return p
            return None

        bin_src = _first_existing(
            find_any(pred_dir_primary, BIN_NAMES_IN),
            find_any(pred_dir_legacy,  BIN_NAMES_IN),
            find_any(predict_vf_dir,   BIN_NAMES_IN),
        )
        multi_src = _first_existing(
            find_any(pred_dir_primary, MULTI_NAMES_IN),
            find_any(pred_dir_legacy,  MULTI_NAMES_IN),
            find_any(predict_vf_dir,   MULTI_NAMES_IN),
        )

        _write_progress(run_dir, 75, "Publishing results…")

        same_dir_for_bin   = bin_src and (os.path.dirname(os.path.abspath(bin_src)) == os.path.abspath(public_pred))
        same_dir_for_multi = multi_src and (os.path.dirname(os.path.abspath(multi_src)) == os.path.abspath(public_pred))
        same_dir_any = bool(same_dir_for_bin or same_dir_for_multi)

        if not same_dir_any:
            clean_dir(public_pred)

        copied = []

        if same_dir_for_bin:
            if os.path.basename(bin_src) != BIN_NAME_OUT:
                os.replace(bin_src, os.path.join(public_pred, BIN_NAME_OUT))
            copied.append("binary")
        elif bin_src:
            copy_with_retry(bin_src, os.path.join(public_pred, BIN_NAME_OUT))
            copied.append("binary")

        if same_dir_for_multi:
            if os.path.basename(multi_src) != MULTI_NAME_OUT:
                os.replace(multi_src, os.path.join(public_pred, MULTI_NAME_OUT))
            copied.append("multiclass")
            try:
                out_multi = os.path.join(public_pred, MULTI_NAME_OUT)
                ids = _ids_from_csv(out_multi, id_col="ID")
                out_fa = os.path.join(public_pred, "multiclass_sequences.fasta")
                _write_subset_fasta(inp_path, ids, out_fa)
            except Exception:
                pass
        elif multi_src:
            out_multi = os.path.join(public_pred, MULTI_NAME_OUT)
            copy_with_retry(multi_src, out_multi)
            copied.append("multiclass")
            try:
                ids = _ids_from_csv(out_multi, id_col="ID")
                out_fa = os.path.join(public_pred, "multiclass_sequences.fasta")
                _write_subset_fasta(inp_path, ids, out_fa)
            except Exception:
                pass

        if os.path.isdir(src_features_dir) and (os.path.abspath(src_features_dir) != os.path.abspath(public_features)):
            for fn in os.listdir(src_features_dir):
                sp = os.path.join(src_features_dir, fn)
                if os.path.isfile(sp):
                    try:
                        copy_with_retry(sp, os.path.join(public_features, fn))
                    except Exception:
                        pass

        if not copied:
            msg = "Completed (no prediction files found)"
        elif "multiclass" not in copied:
            msg = "Completed (no multiclass predictions)"
        else:
            msg = "Completed"

        _write_progress(run_dir, 100, msg, done=True, redirect=reverse("vf:vf_outputs"))

    except subprocess.CalledProcessError as e:
        _write_progress(run_dir, 100, "Pipeline failed", done=True, error=str(e))
    except Exception as e:
        _write_progress(run_dir, 100, "Error", done=True, error=str(e))


# ============================== Views ======================================

def vf_predict_view(request):
    if request.method == "POST":
        form = UploadFastaForm(request.POST, request.FILES)
        if not form.is_valid():
            return render(request, "vf/upload_fasta.html", {"form": form, "error": "❌ Invalid form submission."})

        uploaded = request.FILES["fasta_file"]
        run_id = str(uuid.uuid4())

        vf_uploads = ensure_dir(os.path.join(settings.MEDIA_ROOT, "vf_uploads"))
        inp_path = os.path.join(vf_uploads, f"input_{run_id}.fasta")
        with open(inp_path, "wb+") as f:
            for chunk in uploaded.chunks():
                f.write(chunk)

        run_dir = ensure_dir(os.path.join(settings.MEDIA_ROOT, "vf_runs", run_id))
        _write_progress(run_dir, 5, "Starting…")
        threading.Thread(target=_run_vf_job, args=(inp_path, run_id), daemon=True).start()

        return redirect(reverse("vf:vf_progress", args=[run_id]))

    return render(request, "vf/upload_fasta.html", {"form": UploadFastaForm()})


def vf_progress_view(request, run_id: str):
    return render(request, "vf/progress.html", {"run_id": str(run_id), "results_url": reverse("vf:vf_outputs")})


def vf_progress_status(request, run_id: str):
    run_dir = os.path.join(settings.MEDIA_ROOT, "vf_runs", str(run_id))
    pj = os.path.join(run_dir, "progress.json")
    fallback = {"percent": 0, "message": "Starting…", "done": False, "error": None, "redirect": reverse("vf:vf_outputs")}
    if not os.path.isfile(pj):
        return JsonResponse(fallback)
    try:
        with open(pj, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return JsonResponse(fallback)
    data.setdefault("redirect", reverse("vf:vf_outputs"))
    return JsonResponse(data)


def vf_status_api(request, run_id: str):
    return vf_progress_status(request, run_id)


# ======================= Results & Lanes Data ==============================

def _probe_dirs(logical_folder: str):
    if logical_folder == "pred_vf":
        return [
            os.path.join(settings.MEDIA_ROOT, "pred_vf"),
            os.path.join(settings.MEDIA_ROOT, "vf_uploads", "pred_vf"),
        ]
    elif logical_folder == "features_vf":
        return [
            os.path.join(settings.MEDIA_ROOT, "features_vf"),
            os.path.join(settings.MEDIA_ROOT, "vf_uploads", "features_vf"),
        ]
    else:
        return [os.path.join(settings.MEDIA_ROOT, logical_folder)]


def _collect_items_for(logical_folder: str):
    HIDE_FEATURE_FILES = set()

    newest_by_name = {}
    for base in _probe_dirs(logical_folder):
        if not os.path.isdir(base):
            continue
        for name in os.listdir(base):
            if logical_folder == "features_vf" and name in HIDE_FEATURE_FILES:
                continue
            abs_path = os.path.join(base, name)
            if not os.path.isfile(abs_path):
                continue
            mt = os.path.getmtime(abs_path)
            prev = newest_by_name.get(name)
            if prev is None or mt > prev[1]:
                newest_by_name[name] = (abs_path, mt)

    if logical_folder == "pred_vf":
        allowed = {
            "binary_prediction.csv",
            "multiclass_only_predictions.csv",
            "multiclass_sequences.fasta",
            "prediction_vf_binary.csv",
            "prediction_vf_multiclass.csv",
        }
        newest_by_name = {n: v for n, v in newest_by_name.items() if n in allowed}

    items = []
    for name, (abs_path, _) in newest_by_name.items():
        size_h = _bytes_to_human(os.path.getsize(abs_path))
        rel = os.path.relpath(abs_path, settings.MEDIA_ROOT).replace("\\", "/")
        items.append({
            "folder": logical_folder,
            "name": name,
            "url": settings.MEDIA_URL + rel,
            "size_h": size_h,
        })

    items.sort(key=lambda x: x["name"].lower())
    return items


def _latest_file(candidates):
    exists = [(p, os.path.getmtime(p)) for p in candidates if p and os.path.isfile(p)]
    if not exists:
        return None
    exists.sort(key=lambda x: x[1], reverse=True)
    return exists[0][0]


def _latest_multiclass_csv():
    places = []
    for base in _probe_dirs("pred_vf"):
        if not os.path.isdir(base):
            continue
        places.append(os.path.join(base, "multiclass_only_predictions.csv"))
        places.append(os.path.join(base, "prediction_vf_multiclass.csv"))  # legacy
    return _latest_file(places)


def _latest_binary_csv():
    places = []
    for base in _probe_dirs("pred_vf"):
        if not os.path.isdir(base):
            continue
        places.append(os.path.join(base, "binary_prediction.csv"))
        places.append(os.path.join(base, "prediction_vf_binary.csv"))  # legacy
    return _latest_file(places)


def _build_function_lanes(multiclass_csv, binary_csv):
    """
    Returns:
      lanes: [{"function": str, "ids":[{"id": str, "label": 0|1|None}, ...]}, ...]
      summary: {"non": int, "vir": int}
    """
    label_by_id = {}
    non_cnt = vir_cnt = 0

    if binary_csv and os.path.isfile(binary_csv):
        with open(binary_csv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames:
                cols_lower = [c.lower() for c in reader.fieldnames]
                id_col = None
                lab_col = None
                for i, c in enumerate(cols_lower):
                    if c in ("id", "seq_id", "sequence_id", "name"):
                        id_col = reader.fieldnames[i]
                    if c in ("label", "binary", "pred", "prediction", "virulence", "class", "target"):
                        lab_col = reader.fieldnames[i]
                if id_col is None:
                    id_col = reader.fieldnames[0]
                if lab_col is None and len(reader.fieldnames) > 1:
                    lab_col = reader.fieldnames[1]

                for row in reader:
                    rid = (row.get(id_col) or "").strip().split()[0]
                    val = (row.get(lab_col) or "").strip()
                    if not rid:
                        continue
                    lab = None
                    if val in ("1", "True", "true", "virulent"):
                        lab = 1
                    elif val in ("0", "False", "false", "non-virulent", "non_virulent"):
                        lab = 0
                    if lab is not None:
                        label_by_id[rid] = lab
                        if lab == 1:
                            vir_cnt += 1
                        else:
                            non_cnt += 1

    lanes_map = {}
    if multiclass_csv and os.path.isfile(multiclass_csv):
        with open(multiclass_csv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames:
                cols_lower = [c.lower() for c in reader.fieldnames]
                id_col = func_col = None
                for i, c in enumerate(cols_lower):
                    if c in ("id", "seq_id", "sequence_id", "name"):
                        id_col = reader.fieldnames[i]
                    if ("class" in c) or ("function" in c) or ("vf" in c):
                        func_col = reader.fieldnames[i]
                if id_col is None:
                    id_col = reader.fieldnames[0]
                if func_col is None:
                    func_col = reader.fieldnames[-1]

                for row in reader:
                    rid = (row.get(id_col) or "").strip().split()[0]
                    func = (row.get(func_col) or "").strip() or "Unassigned"
                    if not rid:
                        continue
                    lab = label_by_id.get(rid, None)
                    lanes_map.setdefault(func, []).append({"id": rid, "label": lab})

    lanes = [{"function": k, "ids": v} for k, v in lanes_map.items()]
    lanes.sort(key=lambda x: (-len(x["ids"]), x["function"].lower()))
    return lanes, {"non": non_cnt, "vir": vir_cnt}


def vf_outputs_view(request):
    pred_items = _collect_items_for("pred_vf")
    feat_items = _collect_items_for("features_vf")

    multi_csv = _latest_multiclass_csv()
    bin_csv = _latest_binary_csv()

    lanes, summary = _build_function_lanes(multi_csv, bin_csv) if (multi_csv or bin_csv) else ([], {"non": 0, "vir": 0})

    zip_links = {
        "features_vf": reverse("vf:vf_outputs_zip", args=["features_vf"]),
        "pred_vf":     reverse("vf:vf_outputs_zip", args=["pred_vf"]),
    }

    return render(request, "vf/results_list.html", {
        "items": pred_items,
        "pred_items": pred_items,
        "feat_items": feat_items,
        "zip_links": zip_links,
        "show_sizes": False,
        "show_features": False,
        "lanes_json": json.dumps(lanes),
        "binary_summary_json": json.dumps(summary),
    })


def vf_outputs_zip(request, folder: str):
    if folder not in ("features_vf", "pred_vf"):
        raise Http404("Unknown folder.")

    abs_files = []
    for base in _probe_dirs(folder):
        if not os.path.isdir(base):
            continue
        for root, _, files in os.walk(base):
            for fn in files:
                if folder == "pred_vf" and fn not in {
                    "binary_prediction.csv",
                    "multiclass_only_predictions.csv",
                    "multiclass_sequences.fasta",
                    "prediction_vf_binary.csv",
                    "prediction_vf_multiclass.csv",
                }:
                    continue
                abs_files.append(os.path.join(root, fn))

    if not abs_files:
        raise Http404("Nothing to zip.")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for fpath in abs_files:
            rel = os.path.relpath(fpath, settings.MEDIA_ROOT).replace("\\", "/")
            z.write(fpath, rel)
    buf.seek(0)
    return FileResponse(buf, as_attachment=True, filename=f"{folder}.zip")
