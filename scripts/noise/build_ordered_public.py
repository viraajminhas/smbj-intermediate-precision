"""Recording-ordered trace caches for the public deposits (numeric file order), same window and tunnelling
filter as the pre-registered analyses (1e-5 < G/G0 < 0.5, >= 8 samples, R^2 < 0.95 dropped). Writes one pickle
per deposit: {session: {"order": [file index], "traces": [log10 G arrays], "kept": [bool per file]}}.
No statistic is computed here. Usage: python build_ordered_public.py nbd|liv|pah"""
import collections
import pickle
import re
import sys
import zipfile

import numpy as np

FLOOR, UPPER, MIN_PTS = 1e-5, 0.5, 8
OUT = r"C:\Users\viraa\AppData\Local\Temp\claude\C--Users-viraa-OneDrive-Documents-Claude-Projects-Research\01bcb72f-cee6-4e15-b0db-32772b9f52a1\scratchpad"


def rows_of(txt):
    out = []
    for line in txt.splitlines():
        f = line.split()
        if len(f) < 2:
            continue
        try:
            out.append((float(f[0]), float(f[1])))
        except ValueError:
            continue
    return np.array(out) if len(out) >= MIN_PTS else None


def screen(x, g):
    """Returns (kept, log10 G samples in window) ; kept False if too short or tunnelling-like (R^2 >= 0.95)."""
    m = (g > FLOOR) & (g < UPPER) & np.isfinite(g)
    xs, y = x[m], np.log10(g[m])
    if len(y) < MIN_PTS:
        return False, None
    A = np.vstack([xs, np.ones_like(xs)]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    ss_res = float(np.sum((y - A @ coef) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return r2 < 0.95, y.astype(np.float32)


def num_key(name):
    """(date, time, N) from scan_YYMMDD_HHMM_N.dat or scanYYMMDD_k_N.dat; falls back to the trailing integer."""
    base = name.split("/")[-1]
    m = re.match(r"scan_?(\d{6})_?(\d*)_(\d+)\.dat$", base)
    if m:
        return (m.group(1), m.group(2), int(m.group(3)))
    m = re.search(r"(\d+)\.(dat|txt)$", base)
    return ("", "", int(m.group(1)) if m else 0)


def delft_sections(txt):
    head, _, rest = txt.partition("@")
    parts = rest.split("@")
    return head, parts


def run_nbd():
    z = zipfile.ZipFile(r"C:\Projects\MCBJ_Delft\Data.zip")
    out = {}
    import datetime
    runs = sorted({"/".join(n.split("/")[:-1]) + "/" for n in z.namelist() if n.endswith(".dat")})
    for run in runs:
        names = sorted([n for n in z.namelist() if n.startswith(run) and n.endswith(".dat")], key=num_key)
        rec = dict(order=[], traces=[], kept=[], files=len(names), mtime=[], file_index=[])
        for n in names:
            head, parts = delft_sections(z.read(n).decode("latin1"))
            a = rows_of(parts[0]) if parts else None      # section 1 = breaking trace, as delft_nbd.py
            if a is None:
                continue
            k, y = screen(a[:, 0], a[:, 1])
            rec["kept"].append(k); rec["file_index"].append(num_key(n)[2])
            if k:
                rec["order"].append(num_key(n)[2]); rec["traces"].append(y)
                rec["mtime"].append(datetime.datetime(*z.getinfo(n).date_time).timestamp())
        out[run] = rec
        print(run, "files", len(names), "kept", len(rec["traces"]), "first/last index", rec["order"][:1], rec["order"][-1:], flush=True)
    pickle.dump(out, open(OUT + r"\ordered_nbd.pkl", "wb"))


def run_liv():
    z = zipfile.ZipFile(r"C:\Projects\STMBJ_Liverpool\MCO.zip")
    out = {}
    for p in ("MnPARA/300 mV/", "MnPARA/300 mV_2/"):
        names = sorted([n for n in z.namelist() if n.startswith(p) and n.split("/")[-1].startswith("combined") and n.endswith(".txt")], key=num_key)
        rec = dict(order=[], traces=[], kept=[], files=len(names))
        for n in names:
            a = rows_of(z.read(n).decode("latin1"))
            if a is None:
                continue
            k, y = screen(a[:, 0], a[:, 1])
            rec["kept"].append(k)
            if k:
                rec["order"].append(num_key(n)[2]); rec["traces"].append(y)
        out[p] = rec
        print(p, "files", len(names), "kept", len(rec["traces"]), flush=True)
    pickle.dump(out, open(OUT + r"\ordered_liv.pkl", "wb"))


def run_pah(cap=10000):
    z = zipfile.ZipFile(r"C:\Projects\MCBJ_Delft\PAH\PAH_RT_MCBJ.zip")
    index = collections.defaultdict(list)
    for n in z.namelist():
        p = n.split("/")
        if n.startswith("PAH_RT_MCBJ/Data/") and n.lower().endswith(".dat") and len(p) >= 7:
            index[(p[2], p[3], p[4], p[5].lower().split(" ")[0])].append(n)
    SESS = [("Piezo", "2-OS", "23-2-2023"), ("Piezo", "2-OS", "sample1"), ("Piezo", "2-OS", "sample2"),
            ("Piezo", "SMe-2OS", "Sample1"), ("Piezo", "SMe-2OS", "sample2"),
            ("LC", "2-OS in-situ", "8-2-2023"), ("LC", "2-OS in-situ", "15-2-2023"), ("LC", "2-OS in-situ", "20-2-2023")]
    out = {}
    for key in SESS:
        files = [f for (s, m, sa, st), fs in index.items() if (s, m, sa) == key and st.startswith("mol") for f in fs]
        files = sorted(files, key=lambda n: (n.rsplit("/", 1)[0],) + num_key(n))   # subfolder, then date, time, N
        rec = dict(order=[], traces=[], kept=[], files=len(files), subfolder=[], bias_skipped=0)
        read = 0
        for n in files:
            if read >= cap:
                break
            head, parts = delft_sections(z.read(n).decode("latin1"))
            sv = re.search(r"set voltage : ([0-9.]+)", head)
            if sv is None or abs(float(sv.group(1)) - 0.10) > 0.005:
                rec["bias_skipped"] += 1
                continue
            a = rows_of("\n".join(parts)) if parts else None     # both sections, as pah_mcbj.py
            if a is None:
                continue
            read += 1
            k, y = screen(a[:, 0], a[:, 1])
            rec["kept"].append(k)
            if k:
                rec["order"].append(read); rec["traces"].append(y); rec["subfolder"].append(n.rsplit("/", 1)[0])
        out["/".join(key)] = rec
        print("/".join(key), "files", len(files), "read at 0.10 V", read, "kept", len(rec["traces"]),
              "subfolders", len(set(rec["subfolder"])), "bias skipped", rec["bias_skipped"], flush=True)
    pickle.dump(out, open(OUT + r"\ordered_pah.pkl", "wb"))


if __name__ == "__main__":
    {"nbd": run_nbd, "liv": run_liv, "pah": run_pah}[sys.argv[1]]()
