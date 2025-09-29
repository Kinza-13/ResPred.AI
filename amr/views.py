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


# ================= IDs from MULTICLASS CSV + FASTA export ==================

_ID_CANDIDATES = ("id", "seq_id", "sequence_id", "name", "header", "accession", "sequence id", "seq id")

def _guess_id_col(fieldnames):
    if not fieldnames:
        return None
    lower = [c.lower().strip() for c in fieldnames]
    for cand in _ID_CANDIDATES:
        for i, h in enumerate(lower):
            if cand == h:
                return fieldnames[i]
    # fallbacks
    for i, h in enumerate(lower):
        if "id" == h or h.endswith("_id"):
            return fieldnames[i]
    return fieldnames[0]


# --- ID normalization utilities --------------------------------------------

_PREFIXES = ("sp|", "tr|", "ref|", "gb|", "emb|", "dbj|", "pir|", "pdb|", "lcl|")

def _norm_one(token: str) -> str:
    t = (token or "").strip()
    if not t:
        return ""
    # remove leading '>' if present
    if t.startswith(">"):
        t = t[1:].strip()
    # keep only the first whitespace-separated token
    t = t.split()[0]
    # lowercase compare
    t = t.lower()
    # drop common prefixes (first piece)
    for p in _PREFIXES:
        if t.startswith(p):
            t = t[len(p):]
            break
    return t

def _variants_from_token(token: str):
    """
    Given a header/CSV token, return a set of plausible variants for matching.
    e.g. 'sp|Q9ABC1|NAME' -> {'sp|q9abc1|name', 'q9abc1|name', 'q9abc1'}
    """
    base = _norm_one(token)
    if not base:
        return set()
    out = {base}
    # split on '|' and keep first and last meaningful pieces
    parts = [p for p in base.split("|") if p]
    if parts:
        out.add(parts[0])
        out.add(parts[-1])
    return out

def _expand_keep_ids(raw_ids: set[str]) -> set[str]:
    expanded = set()
    for rid in raw_ids:
        expanded |= _variants_from_token(rid)
    return expanded


def _ids_from_multiclass_csv(csv_path: str) -> set[str]:
    ids = set()
    if not os.path.isfile(csv_path):
        return ids
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return ids
        id_col = _guess_id_col(reader.fieldnames)
        for row in reader:
            rid = (row.get(id_col) or "").strip()
            if not rid:
                continue
            # keep the raw token; we’ll expand into variants for matching
            ids.add(rid)
    return ids


def _write_subset_fasta(src_fasta: str, keep_ids_raw: set[str], out_fasta: str, log_path: str | None = None) -> int:
    """
    Export only sequences whose header ID matches any normalized CSV ID variant.
    Compares case-insensitively and tolerates common prefixes and pipes.
    """
    if not keep_ids_raw or not os.path.isfile(src_fasta):
        # optional log
        if log_path:
            try:
                with open(log_path, "w", encoding="utf-8") as lf:
                    lf.write("No keep IDs or source FASTA missing.\n")
            except Exception:
                pass
        return 0

    keep_ids = _expand_keep_ids(keep_ids_raw)
    written = 0
    total_headers = 0
    matched_headers = 0

    with open(src_fasta, "r", encoding="utf-8") as fi, open(out_fasta, "w", encoding="utf-8") as fo:
        cur_header = None
        cur_id_variants = set()
        cur_seq = []

        def flush():
            nonlocal written, matched_headers
            if cur_header and (cur_id_variants & keep_ids):
                fo.write(cur_header)
                fo.writelines(cur_seq)
                if not cur_seq or not cur_seq[-1].endswith("\n"):
                    fo.write("\n")
                written += 1
                matched_headers += 1

        for line in fi:
            if line.startswith(">"):
                if cur_header is not None:
                    flush()
                cur_header = line
                total_headers += 1
                token = line[1:].strip().split()[0] if len(line) > 1 else ""
                cur_id_variants = _variants_from_token(token)
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_header is not None:
            flush()

    # Remove empty file to avoid a misleading row in the UI
    if written == 0 and os.path.isfile(out_fasta):
        try:
            os.remove(out_fasta)
        except Exception:
            pass

    # optional log for troubleshooting
    if log_path:
        try:
            with open(log_path, "w", encoding="utf-8") as lf:
                lf.write(
                    f"FASTA export report\n"
                    f"  Total FASTA headers: {total_headers}\n"
                    f"  CSV raw IDs: {len(keep_ids_raw)}\n"
                    f"  Expanded ID variants: {len(keep_ids)}\n"
                    f"  Matched headers: {matched_headers}\n"
                    f"  Sequences written: {written}\n"
                )
        except Exception:
            pass

    return written


# ========================== Background job =================================

def _run_amr_job(inp_path: str, run_id: str):
    """
    Run AMR_pred.py; publish outputs to media/features_amr & media/predict_amr.
    Export ONLY: amr_multiclass_sequences.fasta from the MULTICLASS CSV IDs.
    """
    run_dir = ensure_dir(os.path.join(settings.MEDIA_ROOT, "amr_runs", run_id))
    public_features = ensure_dir(os.path.join(settings.MEDIA_ROOT, "features_amr"))
    public_pred     = ensure_dir(os.path.join(settings.MEDIA_ROOT, "predict_amr"))

    try:
        _write_progress(run_dir, 10, "Queued…"); time.sleep(0.15)
        _write_progress(run_dir, 25, "Extracting features…")
        _write_progress(run_dir, 45, "Computing embeddings…")
        _write_progress(run_dir, 60, "Running AMR models…")

        script_path = os.path.join(settings.BASE_DIR, "amr", "AMR_pred.py")

        env = os.environ.copy()
        env["RUN_ID"] = str(run_id)
        env["MEDIA_ROOT"] = str(settings.MEDIA_ROOT)

        subprocess.run([sys.executable, script_path, inp_path],
                       check=True, cwd=settings.BASE_DIR, env=env)

        # Where files may appear
        pref_features_dir   = os.path.join(settings.MEDIA_ROOT, "features_amr")
        legacy_features_dir = os.path.join(settings.MEDIA_ROOT, "amr_uploads", "features_amr")
        src_features_dir = pref_features_dir if (os.path.isdir(pref_features_dir) and os.listdir(pref_features_dir)) else legacy_features_dir

        base_pred_dir   = os.path.join(settings.MEDIA_ROOT, "predict_amr")
        alt_pred_dir    = os.path.join(settings.MEDIA_ROOT, "pred_amr")
        legacy_pred_dir = os.path.join(settings.MEDIA_ROOT, "amr_uploads", "pred_amr")

        bin_name   = "prediction_amr_binary.csv"
        multi_name = "prediction_amr_multiclass.csv"

        bin_src = _first_existing(
            os.path.join(base_pred_dir, bin_name),
            os.path.join(alt_pred_dir,  bin_name),
            os.path.join(legacy_pred_dir, bin_name),
        )
        multi_src = _first_existing(
            os.path.join(base_pred_dir, multi_name),
            os.path.join(alt_pred_dir,  multi_name),
            os.path.join(legacy_pred_dir, multi_name),
        )

        _write_progress(run_dir, 82, "Publishing results…")

        # Publish predictions
        clean_dir(public_pred)
        if bin_src:
            copy_with_retry(bin_src, os.path.join(public_pred, bin_name))
        if multi_src:
            copy_with_retry(multi_src, os.path.join(public_pred, multi_name))

        # ---- FASTA export from MULTICLASS ONLY ----
        try:
            if multi_src:
                ids_multi_raw = _ids_from_multiclass_csv(multi_src)
                log_path = os.path.join(public_pred, "amr_multiclass_fasta_export.log")
                _write_subset_fasta(
                    inp_path,
                    ids_multi_raw,
                    os.path.join(public_pred, "amr_multiclass_sequences.fasta"),
                    log_path=log_path,
                )
        except Exception:
            # do not fail run if export fails
            pass

        # Copy features (best-effort)
        if os.path.isdir(src_features_dir):
            for fn in os.listdir(src_features_dir):
                sp = os.path.join(src_features_dir, fn)
                if os.path.isfile(sp):
                    try:
                        copy_with_retry(sp, os.path.join(public_features, fn))
                    except Exception:
                        pass

        _write_progress(run_dir, 100, "Completed", done=True, redirect=reverse("amr:amr_outputs"))

    except subprocess.CalledProcessError as e:
        _write_progress(run_dir, 100, "Pipeline failed", done=True, error=str(e))
    except Exception as e:
        _write_progress(run_dir, 100, "Error", done=True, error=str(e))


# ============================== Views ======================================

def amr_predict_view(request):
    if request.method == "POST":
        form = UploadFastaForm(request.POST, request.FILES)
        if not form.is_valid():
            return render(request, "amr/upload_fasta.html", {"form": form, "error": "❌ Invalid form submission."})

        uploaded = request.FILES["fasta_file"]
        run_id = str(uuid.uuid4())

        amr_uploads = ensure_dir(os.path.join(settings.MEDIA_ROOT, "amr_uploads"))
        inp_path = os.path.join(amr_uploads, f"input_{run_id}.fasta")
        with open(inp_path, "wb+") as f:
            for chunk in uploaded.chunks():
                f.write(chunk)

        run_dir = ensure_dir(os.path.join(settings.MEDIA_ROOT, "amr_runs", run_id))
        _write_progress(run_dir, 5, "Starting…")
        threading.Thread(target=_run_amr_job, args=(inp_path, run_id), daemon=True).start()

        return redirect(reverse("amr:amr_progress", args=[run_id]))

    return render(request, "amr/upload_fasta.html", {"form": UploadFastaForm()})


def amr_progress_view(request, run_id: str):
    return render(request, "amr/progress.html", {"run_id": str(run_id), "results_url": reverse("amr:amr_outputs")})


def amr_progress_status(request, run_id: str):
    run_dir = os.path.join(settings.MEDIA_ROOT, "amr_runs", str(run_id))
    pj = os.path.join(run_dir, "progress.json")
    fallback = {"percent": 0, "message": "Starting…", "done": False, "error": None, "redirect": reverse("amr:amr_outputs")}
    if not os.path.isfile(pj):
        return JsonResponse(fallback)
    try:
        with open(pj, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return JsonResponse(fallback)
    data.setdefault("redirect", reverse("amr:amr_outputs"))
    return JsonResponse(data)


def amr_status_api(request, run_id: str):
    return amr_progress_status(request, run_id)


# ======================= Results & Lanes Data ==============================

def _probe_dirs(logical_folder: str):
    if logical_folder == "pred_amr":
        return [
            os.path.join(settings.MEDIA_ROOT, "predict_amr"),
            os.path.join(settings.MEDIA_ROOT, "pred_amr"),
            os.path.join(settings.MEDIA_ROOT, "amr_uploads", "pred_amr"),
        ]
    elif logical_folder == "features_amr":
        return [
            os.path.join(settings.MEDIA_ROOT, "features_amr"),
            os.path.join(settings.MEDIA_ROOT, "amr_uploads", "features_amr"),
        ]
    else:
        return [os.path.join(settings.MEDIA_ROOT, logical_folder)]


def _collect_items_for(logical_folder: str):
    HIDE_FEATURE_FILES = {"selected_feat_amr.csv"} if logical_folder == "features_amr" else set()

    newest_by_name = {}
    for base in _probe_dirs(logical_folder):
        if not os.path.isdir(base):
            continue
        for name in os.listdir(base):
            if name in HIDE_FEATURE_FILES:
                continue
            abs_path = os.path.join(base, name)
            if not os.path.isfile(abs_path):
                continue
            mt = os.path.getmtime(abs_path)
            prev = newest_by_name.get(name)
            if prev is None or mt > prev[1]:
                newest_by_name[name] = (abs_path, mt)

    if logical_folder == "pred_amr":
        allowed = {
            "prediction_amr_binary.csv",
            "prediction_amr_multiclass.csv",
            "amr_multiclass_sequences.fasta",
            "amr_multiclass_fasta_export.log",  # expose the tiny log for debugging (optional)
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


def _latest_file(paths):
    candidates = [(p, os.path.getmtime(p)) for p in paths if p and os.path.isfile(p)]
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[1], reverse=True)
    return candidates[0][0]


def _latest_multiclass_csv():
    paths = []
    for base in _probe_dirs("pred_amr"):
        if not os.path.isdir(base):
            continue
        paths.append(os.path.join(base, "prediction_amr_multiclass.csv"))
    return _latest_file(paths)


def _latest_binary_csv():
    paths = []
    for base in _probe_dirs("pred_amr"):
        if not os.path.isdir(base):
            continue
        paths.append(os.path.join(base, "prediction_amr_binary.csv"))
    return _latest_file(paths)


def _build_family_lanes(multiclass_csv: str | None, binary_csv: str | None):
    """
    Lanes for the UI (kept as before).
    Uses binary for colouring if present; otherwise shows unlabelled tiles.
    """
    id2lab = {}
    non_cnt = vir_cnt = 0

    if binary_csv and os.path.isfile(binary_csv):
        with open(binary_csv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames:
                # id column
                lower = [c.lower() for c in reader.fieldnames]
                id_col = None; lab_col = None
                for i, h in enumerate(lower):
                    if h in ("id", "seq_id", "sequence_id", "name", "header", "accession", "sequence id", "seq id"):
                        id_col = reader.fieldnames[i]
                    if any(k in h for k in ("label","binary","class","resist","pred","prediction","target")):
                        lab_col = reader.fieldnames[i]
                if id_col is None: id_col = reader.fieldnames[0]
                if lab_col is None and len(reader.fieldnames) > 1: lab_col = reader.fieldnames[1]
                for row in reader:
                    rid = (row.get(id_col) or "").strip().split()[0]
                    val = (row.get(lab_col) or "").strip().lower()
                    if not rid:
                        continue
                    if val in ("1","true","resistant","yes","virulent"): id2lab[rid]=1; vir_cnt+=1
                    elif val in ("0","false","non-resistant","susceptible","no","non_virulent"): id2lab[rid]=0; non_cnt+=1

    lanes_map = {}
    if multiclass_csv and os.path.isfile(multiclass_csv):
        with open(multiclass_csv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames:
                id_col = _guess_id_col(reader.fieldnames)
                lower = [c.lower() for c in reader.fieldnames]
                fam_col = None
                for i, h in enumerate(lower):
                    if any(k in h for k in ("family","class","amr","gene_family","label")):
                        fam_col = reader.fieldnames[i]; break
                if fam_col is None: fam_col = reader.fieldnames[-1]
                for row in reader:
                    rid = (row.get(id_col) or "").strip().split()[0]
                    fam = (row.get(fam_col) or "").strip() or "Unassigned"
                    if not rid:
                        continue
                    lab = id2lab.get(rid, None)
                    lanes_map.setdefault(fam, []).append({"id": rid, "label": lab})

    lanes = [{"function": k, "ids": v} for k, v in lanes_map.items()]
    lanes.sort(key=lambda x: (-len(x["ids"]), x["function"].lower()))
    return lanes, {"non": non_cnt, "vir": vir_cnt}


def amr_outputs_view(request):
    items = _collect_items_for("pred_amr")

    multi_csv = _latest_multiclass_csv()
    bin_csv   = _latest_binary_csv()
    lanes, summary = _build_family_lanes(multi_csv, bin_csv) if (multi_csv or bin_csv) else ([], {"non": 0, "vir": 0})

    zip_links = {
        "features_amr": reverse("amr:amr_outputs_zip", args=["features_amr"]),
        "pred_amr":     reverse("amr:amr_outputs_zip", args=["pred_amr"]),
    }

    return render(request, "amr/result_lists.html", {
        "items": items,
        "zip_links": zip_links,
        "show_sizes": False,
        "show_features": False,
        "lanes_json": json.dumps(lanes),
        "binary_summary_json": json.dumps(summary),
    })


def amr_outputs_zip(request, folder: str):
    if folder not in ("features_amr", "pred_amr"):
        raise Http404("Unknown folder.")

    abs_files = []
    for base in _probe_dirs(folder):
        if not os.path.isdir(base):
            continue
        for root, _, files in os.walk(base):
            for fn in files:
                if folder == "pred_amr" and fn not in {
                    "prediction_amr_binary.csv",
                    "prediction_amr_multiclass.csv",
                    "amr_multiclass_sequences.fasta",
                    "amr_multiclass_fasta_export.log",
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
