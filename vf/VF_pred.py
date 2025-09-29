import os
import sys
import json
import csv
import re

import pandas as pd
import requests
from Bio import SeqIO
import peptides

# ==== Live progress bridge (writes to media/amr_runs/<RUN_ID>/progress.json) ====
RUN_ID = os.environ.get("RUN_ID")
MEDIA_ROOT = os.environ.get("MEDIA_ROOT", "media")
RUN_DIR = os.path.join(MEDIA_ROOT, "vf_runs", RUN_ID) if RUN_ID else None
PROG_PATH = os.path.join(RUN_DIR, "progress.json") if RUN_DIR else None
if RUN_DIR:
    os.makedirs(RUN_DIR, exist_ok=True)

def _progress(pct=None, msg=None, done=None, error=None):
    """Merge-update progress.json (atomic write). Non-fatal if file missing."""
    if not PROG_PATH:
        return
    try:
        data = {"percent": 0, "message": "Starting…", "done": False, "error": None, "redirect": None}
        if os.path.exists(PROG_PATH):
            try:
                with open(PROG_PATH, "r", encoding="utf-8") as f:
                    data.update(json.load(f))
            except Exception:
                pass
        if pct   is not None: data["percent"] = int(pct)
        if msg   is not None: data["message"] = str(msg)
        if done  is not None: data["done"] = bool(done)
        if error is not None: data["error"] = str(error)
        tmp = PROG_PATH + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(data, f)
        os.replace(tmp, PROG_PATH)
    except Exception:
        pass
# ==== /progress bridge ====

# === Output directories ===
# Write to a staging area (amr_uploads/*) so views.py can publish/copy into public folders.
WORK_BASE    = os.path.join(MEDIA_ROOT, "vf_uploads")
FEATURES_DIR = os.path.join(WORK_BASE, "features_vf")
PRED_DIR     = os.path.join(WORK_BASE, "pred_vf")
os.makedirs(FEATURES_DIR, exist_ok=True)
os.makedirs(PRED_DIR, exist_ok=True)

# === Inputs / filenames ===
fasta_file       = sys.argv[1] if len(sys.argv) > 1 else "deploy_amr.fasta"
peptide_csv      = os.path.join(FEATURES_DIR, "physicochemical_prop_vf.csv")
ptm_raw_csv      = os.path.join(FEATURES_DIR, "post_translational_modifi_vf.csv")
ptm_avg_csv      = os.path.join(FEATURES_DIR, "average_ptms_vf.csv")
final_output_csv = os.path.join(FEATURES_DIR, "extract_feat_vf.csv")
esm_output_tsv   = os.path.join(FEATURES_DIR, "combinefea_12_esm2.tsv")  # aligned name
cutoff = 0.5

# total sequences for progress increments
try:
    _total_seqs = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
except Exception:
    _total_seqs = None

ptm_models = [
    "Phosphoserine_Phosphothreonine",
    "Phosphotyrosine",
    "N-linked_glycosylation",
    "O-linked_glycosylation",
    "Methylarginine",
    "N6-acetyllysine",
    "Ubiquitination",
    "Methyllysine",
    "Pyrrolidone_carboxylic_acid",
    "S-palmitoyl_cysteine",
    "Hydroxyproline",
    "Hydroxylysine",
    "SUMOylation",
]
model_string = ";".join(ptm_models)

# === STEP 1: PEPTIDE FEATURE EXTRACTION (matches original schema) ===
## === STEP 1: PEPTIDE FEATURE EXTRACTION (target schema with ID; no AF/PRIN/VSTPV) ===
print(" Extracting peptide features...")
_progress(25, "Extracting Physicochemical features…")

peptide_data_list = []

for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
    seq = str(record.seq)
    pep = peptides.Peptide(seq)

    # Core metrics
    row = {
        "ID":                 record.id,                  # <-- ID (not "Sequence ID")
        "Aliphatic Index":    pep.aliphatic_index(),
        "Boman Index":        pep.boman(),
        "Charge at pH 7.4":   pep.charge(pH=7.4),
        "Isoelectric Point":  pep.isoelectric_point(),
        "Hydrophobicity":     pep.hydrophobicity(),
        "Hydrophobic Moment": pep.hydrophobic_moment(),
        "Instability Index":  pep.instability_index(),
        "Molecular Weight":   pep.molecular_weight(),
    }

    # Vector blocks with the verbose names (e.g., "BLOSUM Indices 1", …)
    blocks = {
        "BLOSUM Indices":      pep.blosum_indices(),
        "Cruciani Properties": pep.cruciani_properties(),
        "FASGAI Vectors":      pep.fasgai_vectors(),
        "Kidera Factors":      pep.kidera_factors(),
        "MS-WHIM Scores":      pep.ms_whim_scores(),   # block form
        "PCP Descriptors":     pep.pcp_descriptors(),
        "ProtFP Descriptors":  pep.protfp_descriptors(),
        "Sneath Vectors":      pep.sneath_vectors(),
        "ST Scales":           pep.st_scales(),
        "T Scales":            pep.t_scales(),
        "VHSE Scales":         pep.vhse_scales(),
        "Z Scales":            pep.z_scales(),
    }
    for name, vec in blocks.items():
        for j, val in enumerate(vec):
            row[f"{name} {j+1}"] = val

    # Individual WHIM fields ("WHIM Score 1..3")
    for j, s in enumerate(pep.ms_whim_scores()):
        row[f"WHIM Score {j+1}"] = s

    # Flat descriptor dict (keeps BLOSUM1.., PP1.., F1.., MSWHIM1.. etc.)
    # Exclude ONLY AF*, PRIN*, VSTPV*
    extra = pep.descriptors()
    for k, v in extra.items():
        k_str = str(k)
        if k_str.startswith(("AF", "PRIN", "VSTPV")):
            continue
        row[k_str] = v

    peptide_data_list.append(row)

    if _total_seqs:
        pct = 25 + int(10 * (i + 1) / _total_seqs)
        if (i + 1) % 2 == 0 or (i + 1) == _total_seqs:
            _progress(pct, f"Physicochemical: {i+1}/{_total_seqs}")

peptide_df = pd.DataFrame(peptide_data_list)
peptide_df.to_csv(peptide_csv, index=False)
print(f" Peptide features saved → {peptide_csv}")
print(f" 🔢 Column count (incl. 'ID'): {peptide_df.shape[1]}")
_progress(35, f"Physicochemical saved → {peptide_csv}")



# === STEP 2: PTM EXTRACTION USING MUSITEDEEP (windowed GET, no skips) ===
print("🔬 Predicting PTMs from MusiteDeep API...")
_progress(50, "Predicting PTMs from MusiteDeep…")

import re, csv, json, time, random
from collections import defaultdict
from urllib.parse import quote
import requests
from Bio import SeqIO

# ---------- CONFIG FOR WINDOWING ----------
WIN = 800        # window length (aa) — increase for fewer API calls
OVERLAP = 32     # overlap (aa) — keeps context across windows
BASE = "https://api.musite.net/musitedeep"

def filter_ptmscore(ptms, cutoff=0.5):
    """Return 'PTM:score' items >= cutoff; 'None' if none."""
    if not ptms or str(ptms).strip() == "":
        return "None"
    out = []
    for item in str(ptms).split(";"):
        if ":" in item:
            k, v = item.split(":", 1)
            try:
                if float(v) >= cutoff:
                    out.append(f"{k}:{v}")
            except ValueError:
                pass
    return ";".join(out) if out else "None"

def _clean_sequence(seq: str) -> str:
    """Uppercase, remove whitespace, replace non-AA with 'X' (tolerated by MusiteDeep)."""
    s = re.sub(r"\s+", "", str(seq).upper())
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", s)

def _windows(n: int, win: int = WIN, overlap: int = OVERLAP):
    """Yield (start,end) 0-based half-open windows covering length n."""
    if n <= win:
        yield (0, n)
        return
    step = win - overlap
    i = 0
    while i < n:
        j = min(i + win, n)
        yield (i, j)
        if j == n:
            break
        i += step

def _parse_ptm_scores(s: str) -> dict:
    """'A:0.5;B:0.2' -> {'A':0.5,'B':0.2} (keep max when duplicated)."""
    d = {}
    for item in str(s).split(";"):
        if ":" in item:
            k, v = item.split(":", 1)
            try:
                v = float(v)
            except ValueError:
                continue
            if v > d.get(k, 0.0):
                d[k] = v
    return d

def _scores_to_string(d: dict) -> str:
    return ";".join(f"{k}:{d[k]}" for k in sorted(d.keys()))

def _musitedeep_get(models: str, seq: str, timeout=25, retries=3):
    """GET request with URL-encoding; tolerant JSON parsing; retries on transient errors."""
    url = f"{BASE}/{quote(models, safe='')}/{quote(seq, safe='')}"
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=timeout)
            if r.ok:
                try:
                    data = r.json()
                except Exception:
                    data = json.loads(r.content.decode("utf-8", errors="ignore"))
                if isinstance(data, dict):
                    return data.get("Results", [])
                if isinstance(data, list):
                    return data
                return []
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(min(2**attempt, 6) + random.uniform(0, 0.5))
                continue
            print(f"❌ MusiteDeep HTTP {r.status_code}")
            return []
        except (requests.Timeout, requests.ConnectionError):
            time.sleep(min(2**attempt, 6) + random.uniform(0, 0.5))
            continue
        except Exception as e:
            print(f"❌ MusiteDeep error: {e}")
            return []
    return []

def _predict_ptms_for_sequence(full_seq: str, models_str: str):
    """
    Split any-length sequence into overlapping windows, call API per window,
    map window positions to global (1-based), and merge overlaps by max score.
    Returns list of dicts like API rows: {'Position','Residue','PTMscores','Cutoff=0.5'}
    """
    seq = _clean_sequence(full_seq)
    n = len(seq)

    per_site = {}  # pos (1-based) -> {'Residue': aa, 'scores': {ptm:score}}

    for (start, end) in _windows(n):
        sub = seq[start:end]
        results = _musitedeep_get(model_string, sub)
        if not isinstance(results, list):
            continue

        for res in results:
            if isinstance(res, dict):
                pos_str = res.get("Position", "")
                residue = res.get("Residue", "")
                ptm_raw = res.get("PTMscores", "")
            else:
                # very rare alternative shapes (be defensive)
                try:
                    pos_str, residue, ptm_raw = str(res[0]), str(res[1]), str(res[2])
                except Exception:
                    continue

            try:
                local_pos = int(str(pos_str).strip())
            except Exception:
                continue

            # API positions are 1-based relative to the sub-sequence
            global_pos = start + local_pos  # keep 1-based indexing

            # fallback residue from full sequence if API didn't send it
            aa = residue or (seq[global_pos - 1] if 1 <= global_pos <= n else "")

            scores = _parse_ptm_scores(ptm_raw)

            if global_pos not in per_site:
                per_site[global_pos] = {"Residue": aa, "scores": {}}
            for k, v in scores.items():
                if v > per_site[global_pos]["scores"].get(k, 0.0):
                    per_site[global_pos]["scores"][k] = v

    # back to a sorted list of rows
    rows = []
    for pos in sorted(per_site.keys()):
        aa = per_site[pos]["Residue"]
        score_str = _scores_to_string(per_site[pos]["scores"])
        rows.append({
            "Position": str(pos),
            "Residue": aa,
            "PTMscores": score_str,
            "Cutoff=0.5": filter_ptmscore(score_str, cutoff)
        })
    return rows

# ---- run over FASTA and write CSV (same format as before) ----
with open(ptm_raw_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["ID", "Position", "Residue", "PTMscores", "Cutoff=0.5"])

    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        seq_id = record.id
        sequence = str(record.seq)
        clean_id = seq_id.strip().lstrip('>')

        # header row per sequence (kept for your averaging step)
        writer.writerow([f">{clean_id}"])

        rows = _predict_ptms_for_sequence(sequence, model_string)

        # write all predicted sites; do not synthesize fake empty rows
        for r in rows:
            writer.writerow([clean_id, r["Position"], r["Residue"], r["PTMscores"], r["Cutoff=0.5"]])

        # progress per sequence
        if _total_seqs:
            pct = 50 + int(20 * (i + 1) / _total_seqs)
            _progress(pct, f"PTMs: {i+1}/{_total_seqs}")

print(f" PTM raw output saved → {ptm_raw_csv}")
_progress(70, f"PTM raw output saved → {ptm_raw_csv}")

# === STEP 3: AVERAGE PTM SCORES (robust) ===
print("Averaging PTM scores...")
_progress(72, "Averaging PTM scores…")

# Read as strings so "None" doesn't become NaN; keep empty strings as empty
df = pd.read_csv(
    ptm_raw_csv,
    dtype=str,
    keep_default_na=False,
    na_filter=False,
)

# Canonical PTM categories (must match MusiteDeep keys)
ptm_categories = [
    "N6-acetyllysine",
    "Methyllysine",
    "SUMOylation",
    "Hydroxylysine",
    "Ubiquitination",
    "N-linked_glycosylation",
    "Phosphothreonine",
    "Phosphoserine",
    "Hydroxyproline",
    "Phosphotyrosine",
    "Pyrrolidone_carboxylic_acid",
    "Methylarginine",
    "S-palmitoyl_cysteine",
    "O-linked_glycosylation",
]

from collections import defaultdict
import math

# sums[id][ptm], counts[id][ptm]
sums   = defaultdict(lambda: {k: 0.0 for k in ptm_categories})
counts = defaultdict(lambda: {k: 0   for k in ptm_categories})

current_id = None

for _, row in df.iterrows():
    rid = str(row.get("ID", "")).strip()

    # Header line marking a new sequence: ">ID"
    if rid.startswith(">"):
        current_id = rid.lstrip(">").strip()
        # Ensure entry exists (even if it ends up with all zeros)
        _ = sums[current_id]; _ = counts[current_id]
        continue

    # Data lines should refer to the same base id
    if not current_id:
        # If file starts without a header (shouldn't), skip safely
        continue

    # Parse PTM scores field safely
    ptm_field = str(row.get("PTMscores", "")).strip()
    if not ptm_field or ptm_field.lower() == "none":
        continue

    # "X:0.123;Y:0.456;..."
    for item in ptm_field.split(";"):
        item = item.strip()
        if not item or ":" not in item:
            continue
        name, val = item.split(":", 1)
        name = name.strip()
        val  = val.strip()
        if name in ptm_categories:
            try:
                x = float(val)
                if not math.isfinite(x):
                    continue
                sums[current_id][name]   += x
                counts[current_id][name] += 1
            except ValueError:
                continue

# Build averages
records = []
for sid in sums.keys():
    rec = {"ID": sid}
    for k in ptm_categories:
        if counts[sid][k] > 0:
            rec[k] = sums[sid][k] / counts[sid][k]
        else:
            rec[k] = 0.0
    records.append(rec)

ptm_avg_df = pd.DataFrame(records)
ptm_avg_df.to_csv(ptm_avg_csv, index=False)

print(f"✅ PTM averages saved → {ptm_avg_csv} (n={len(ptm_avg_df)})")
_progress(74, f"PTM averages saved → {ptm_avg_csv}")

# === STEP 4: ESM2 ENSEMBLE EMBEDDINGS EXTRACTION ===
print("Extracting ESM2 embeddings...")
_progress(75, "Extracting ESM2 embeddings…")  # matches 'ESM2|Embedding'

from transformers import EsmTokenizer, EsmModel
import torch
import numpy as np

class ESM2Embedding:
    def __init__(self, model_name="facebook/esm2_t30_150M_UR50D"):
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmModel.from_pretrained(model_name)
        self.model.eval()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)
        print(f"✅ ESM2 model loaded on {self.device}")

    def _calculate_features(self, sequences):
        embeddings = []
        total = len(sequences)
        for i, seq in enumerate(sequences):
            inputs = self.tokenizer(seq, return_tensors="pt", padding=True, truncation=True, max_length=1024)
            inputs = {key: val.to(self.device) for key, val in inputs.items()}
            with torch.no_grad():
                outputs = self.model(**inputs)
            last_hidden_states = outputs.last_hidden_state
            emb = last_hidden_states.mean(dim=1).detach().cpu().numpy()
            embeddings.append(emb[0])
            if total > 0:
                pct = 75 + int(15 * (i + 1) / total)
                _progress(pct, f"ESM2 embeddings: {i+1}/{total}")
        return np.array(embeddings)

    def process_fasta(self, fasta_file, output_file):
        sequences, ids = self._read_fasta(fasta_file)
        embeddings = self._calculate_features(sequences)
        with open(output_file, "w") as f:
            header = ["ID"] + [f"feature_{i}" for i in range(embeddings.shape[1])]
            f.write("\t".join(header) + "\n")
            for i, emb in enumerate(embeddings):
                f.write(f"{ids[i]}\t" + "\t".join(map(str, emb)) + "\n")
        print(f"✅ ESM2 embeddings saved → {output_file}")

    def _read_fasta(self, fasta_file):
        sequences, ids = [], []
        with open(fasta_file, "r") as f:
            seq_id, seq = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq_id:
                        sequences.append("".join(seq))
                        ids.append(seq_id)
                    seq_id = line[1:].split()[0]
                    seq = []
                else:
                    seq.append(line.upper())
            if seq_id:
                sequences.append("".join(seq))
                ids.append(seq_id)
        return sequences, ids

esm2 = ESM2Embedding()
esm2.process_fasta(fasta_file, esm_output_tsv)
_progress(90, f"ESM2 embeddings saved → {esm_output_tsv}")

# === STEP 5: Merge features ===
print("Merging peptide, PTM average, and ESM2 embeddings (by common ID)...")
_progress(92, "Merging features…")

esm_df = pd.read_csv(esm_output_tsv, sep="\t")
peptide_df = pd.read_csv(peptide_csv)
ptm_avg_df = pd.read_csv(ptm_avg_csv)

esm_df["ID"] = esm_df["ID"].astype(str).str.strip().str.lstrip('>')

if "Sequence ID" in peptide_df.columns and "ID" not in peptide_df.columns:
    peptide_df.rename(columns={"Sequence ID": "ID"}, inplace=True)
peptide_df["ID"] = peptide_df["ID"].astype(str).str.strip().str.lstrip('>')

ptm_avg_df["ID"] = ptm_avg_df["ID"].astype(str).str.strip().str.lstrip('>')

final_df = pd.merge(peptide_df, ptm_avg_df, on="ID", how="inner", suffixes=('_pep', '_ptm'))
final_df = pd.merge(final_df, esm_df, on="ID", how="inner")
final_df.to_csv(final_output_csv, index=False)

print("\n📊 Rows in final merged dataset:", len(final_df))
print(f"✅ All features saved to → {final_output_csv}")



#Feature_selection_vf
import pandas as pd

# === Step 1: Load the dataset
df = pd.read_csv(os.path.join(FEATURES_DIR, "extract_feat_vf.csv"))  # was media/features_vf/...

# === Step 2: Define columns to keep
columns_to_keep = [
    'feature_388','feature_456','feature_373','feature_432','feature_176','feature_242','feature_613','feature_230',
    'feature_473','feature_95','feature_155','feature_203','feature_548','feature_367','feature_542','feature_4',
    'feature_529','feature_28','feature_517','feature_491','feature_199','feature_338','feature_224','feature_592',
    'feature_222','feature_126','feature_525','feature_497','feature_192','feature_506','feature_639','feature_23',
    'feature_236','feature_527','feature_425','feature_289','feature_177','feature_299','feature_80','feature_583',
    'feature_463','feature_99','feature_15','feature_599','feature_437','feature_325','feature_36','feature_302',
    'feature_164','feature_21','feature_217','feature_179','feature_249','feature_290','feature_184','feature_206',
    'feature_52','feature_334','feature_282','feature_364','feature_618','feature_269','feature_391','feature_291',
    'feature_485','feature_154','feature_185','feature_247','feature_554','feature_535','feature_591','feature_20',
    'feature_218','feature_163','feature_227','feature_452','feature_42','feature_423','feature_536','feature_113',
    'feature_85','feature_550','feature_366','feature_311','feature_444','feature_250','feature_454','feature_51',
    'feature_306','feature_305','feature_122','feature_547','feature_568','feature_178','feature_114','feature_275',
    'feature_228','feature_415','feature_226','feature_10','feature_412','feature_386','feature_276','feature_574',
    'feature_87','feature_25','feature_503','feature_598','feature_372','feature_619','feature_221','feature_257',
    'feature_49','feature_213','feature_37','feature_570','feature_30','feature_624','feature_494','feature_354',
    'feature_445','feature_181','feature_251','feature_297','feature_397','feature_194','feature_424','feature_314',
    'feature_460','feature_180','feature_470','feature_158','feature_32','feature_530','feature_170','feature_442',
    'feature_216','feature_34','feature_144','feature_539','feature_523','feature_280','feature_526','feature_337',
    'feature_626','feature_82','feature_380','Isoelectric Point','feature_421','feature_537','feature_458',
    'feature_575','feature_204','feature_351','feature_500','feature_293','feature_369','feature_29','feature_6',
    'feature_630','feature_577','feature_403','feature_435','feature_211','feature_142','feature_519','feature_303',
    'feature_566','feature_109','feature_427','feature_410','feature_515','feature_595','feature_75','feature_121',
    'feature_320','feature_482','SVGER1','feature_576','feature_108','feature_602','feature_471','feature_322',
    'feature_420','feature_316','feature_528','feature_446','feature_119','KF6','feature_438','feature_89',
    'feature_579','feature_492','feature_241','feature_426','feature_182','feature_244','feature_146','feature_324',
    'feature_559','feature_449','feature_159','SVGER7','feature_136','feature_65','feature_551','feature_429',
    'feature_279','SVGER5','feature_587','feature_106','feature_582','feature_104','feature_620','feature_612',
    'feature_631','feature_239','feature_414','feature_2','feature_207','feature_59','Kidera Factors 6','feature_14',
    'feature_26','feature_475','feature_286','feature_356','feature_229','feature_533','feature_215','feature_510',
    'feature_571','feature_57','feature_384','feature_544','feature_22','feature_586','feature_365','feature_349',
    'feature_56','ST Scales 7','feature_16','feature_292','feature_214','feature_370','feature_567','feature_628',
    'feature_231','feature_513','Charge at pH 7.4','feature_260','feature_368','feature_173','feature_309',
    'feature_362','feature_160','feature_498','feature_331','feature_18','feature_198','feature_125','feature_264',
    'MS-WHIM Scores 2','feature_101','MSWHIM2','feature_374','feature_107','feature_259','feature_347','feature_111',
    'Phosphoserine','feature_581','feature_97','feature_219','feature_469','feature_271','feature_404','feature_243',
    'ProtFP Descriptors 7','feature_433','feature_371','feature_27','feature_493','WHIM Score 2','feature_382',
    'feature_235','feature_248','feature_274','feature_105','feature_387','feature_0','feature_150','feature_195',
    'ID'
]

# === Step 3: Filter only selected features
_progress(95, "Selecting features…")
df_selected = df[columns_to_keep]

# === Step 4: Save filtered features to new file
df_selected.to_csv(os.path.join(FEATURES_DIR, "selected_feat_vf.csv"), index=False)
_progress(96, f"Selected features saved → {os.path.join(FEATURES_DIR, 'selected_feat_vf.csv')}")

# === VF_predictions

import pandas as pd
import joblib
import os

# === CONFIGURATION ===
unseen_file = os.path.join(FEATURES_DIR, "selected_feat_vf.csv")
binary_scaler_path = "vf/models_pkl_files_vf/scaler_vf_binary.pkl"
binary_model_path = "vf/models_pkl_files_vf/stacking_model_vf_binary.pkl"

multiclass_scaler_path = "vf/models_pkl_files_vf/scaler_vf_multiclass.pkl"
multiclass_model_path = "vf/models_pkl_files_vf/svc_vf_multiclass.pkl"
label_encoder_path = "vf/models_pkl_files_vf/label_encoder_vf_multiclass.pkl"  # This must be fitted during training


# === STEP 1: Load Unseen Data ===
df_unseen = pd.read_csv(unseen_file)
id_cols = ['ID']  # You may adjust this based on your actual columns
X_unseen = df_unseen.drop(columns=id_cols, errors='ignore')

# === STEP 2: Binary Prediction ===
_progress(97, "Running trained models (binary)…")
binary_scaler = joblib.load(binary_scaler_path)
X_unseen_scaled = binary_scaler.transform(X_unseen)
binary_model = joblib.load(binary_model_path)
binary_preds = binary_model.predict(X_unseen_scaled)

df_unseen["Binary_Prediction"] = binary_preds

# === STEP 3: Save Binary Predictions ===
# === STEP 3: Save Binary Predictions (only ID + Binary_Prediction) ===
binary_out_path = os.path.join(PRED_DIR, "binary_prediction.csv")

# Keep df_unseen unchanged for later steps; just create a slim view for output
id_keep = [c for c in id_cols if c in df_unseen.columns]  # usually ["ID"]
df_binary_out = df_unseen.loc[:, id_keep + ["Binary_Prediction"]]

df_binary_out.to_csv(binary_out_path, index=False)
print(f"✅ Binary prediction file saved to: {binary_out_path}")
_progress(98, f"Binary prediction file saved to: {binary_out_path}")


# === STEP 4: Filter Only Positive (1) Rows ===
df_positive = df_unseen[df_unseen["Binary_Prediction"] == 1].copy()

if df_positive.empty:
    print("⚠️ No positive samples (1s) found. Skipping multiclass prediction.")
    _progress(98, "No positive samples (1s) found. Skipping multiclass prediction.")
else:
    X_positive = df_positive.drop(columns=id_cols + ["Binary_Prediction"], errors='ignore')

    # === STEP 5: Multiclass Prediction ===
    _progress(98, "Running trained models (multiclass)…")
    multiclass_scaler = joblib.load(multiclass_scaler_path)
    X_positive_scaled = multiclass_scaler.transform(X_positive)

    multiclass_model = joblib.load(multiclass_model_path)
    multiclass_preds = multiclass_model.predict(X_positive_scaled)

    # === STEP 6: Decode Categorical Labels ===
    label_encoder = joblib.load(label_encoder_path)
    multiclass_decoded = label_encoder.inverse_transform(multiclass_preds)

    df_positive["Functions_Prediction"] = multiclass_decoded

    # === STEP 7: Save Only Multiclass Rows (Binary = 1) ===
    multiclass_out_path = os.path.join(PRED_DIR, "multiclass_only_predictions.csv")
    df_positive_final = df_positive[id_cols + ["Functions_Prediction"]]
    df_positive_final.to_csv(multiclass_out_path, index=False)
    print(f"✅ Multiclass predictions saved: {multiclass_out_path}")
    _progress(99, f"Multiclass predictions saved: {multiclass_out_path}")
