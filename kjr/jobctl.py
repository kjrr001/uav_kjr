#!/usr/bin/env python3
"""uav_kjr job control: run blocks from control files (*.cfg).

Control file: plain text, one "KEY : value" per line, '#' starts a comment, keys are not case
sensitive. Relative paths are relative to the file they are written in; ${KEY} is replaced by
the value of another key (e.g. WORK_DIR : ../runs/${BLOCK}); INCLUDE : common.cfg reads another
file first (its values are overridden by the including file). Unknown keys are an error.
`kjr_run --keys` lists every key with its default; `kjr_run --stages` the stages.

START .. END (name or number) is the stage range. Stages before START must be complete;
stages in the range are run according to EXISTING:
  skip     a complete stage (marker or its outputs present) is not run again; an incomplete
           one runs and resumes from its partial output; a stale one is cleared and run (default)
  archive  run again; its old outputs are moved to <WORK_DIR>/old/<time>/ first
  delete   run again; its old outputs are deleted first
When a stage runs, every later stage is marked stale (job/stale) until it runs again: a stale
stage is never skipped, and a stale stage before START stops the job.

Every run: <WORK_DIR>/job/<time>/{control.cfg (copy), resolved.cfg (all values), code/ (the repo
as run: editing the repo while a job runs does not change it), git.txt}, per-stage logs in
<WORK_DIR>/job/logs/, stage summary in <WORK_DIR>/kjr_run.log, completion markers job/done/.
"""
import argparse, datetime, difflib, json, os, re, shutil, subprocess, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (name, description) in run order; the stage number is the position (1-based)
STAGES = [
    ("SFM_FEATURE", "OpenSfM: EXIF metadata, camera model (brown from the EXIF focal), SIFT features"),
    ("SFM_MATCH", "OpenSfM: feature matching, tracks"),
    ("SFM_RECON", "OpenSfM: incremental reconstruction, bundle adjustment with the GCPs"),
    ("GEOREF", "OpenSfM ENU -> UTM, ortho grid, reference DSM from the SfM points (osfm_import.py)"),
    ("PAIRS", "stereo pair list (and, REF_DSM : stereo, a reference DSM from the SfM cameras)"),
    ("TIEPOINT", "sub-pixel ties on the original images: SIFT -> MAGSAC++ -> ALSC (GPU)"),
    ("TRACKS", "multi-view tracks re-measured by ALSC (only with USE_TRACKS : yes)"),
    ("BBA", "own bundle block adjustment, radial + 2-D grid image correction"),
    ("UNDISTORT", "Lanczos undistortion with the BBA camera"),
    ("DSM_EXTRACTION", "census-SGM DSM per stereo pair"),
    ("FUSE", "robust fusion of the pair DSMs"),
    ("BLEND", "long waves (> BLEND_SIGMA) from the reference DSM, GCP offset"),
    ("ORTHO", "orthophoto from the final DSM"),
    ("PRODUCTS", "COG / ENVI / JPG deliverables in OUTPUT_DIR"),
    ("DIAG", "colour hillshade, GCP heights, image-locked pattern test"),
    ("EXPORT", "delivery folder EXPORT_DIR/<BLOCK> (small JPGs, footprints, context map); only on request"),
]
NAMES = [s[0] for s in STAGES]

# key: (default, required, is_path, help).  Defaults may use ${KEY}.
KEYS = {
    "TITLE": ("", False, False, "free text, written to the log"),
    "BLOCK": ("", True, False, "block name (product file names, logs)"),
    "INPUT_DIR": ("", True, True, "camera images (*.jpg)"),
    "GCP_FILE": ("", True, True, "GCP list, OpenSfM/ODM format (projection line, then E N Z px py image [id])"),
    "WORK_DIR": ("", True, True, "intermediate results"),
    "OUTPUT_DIR": ("", True, True, "products (DSM, ortho, JPGs)"),
    "START": ("", False, False, "first stage (name or number); default: the first stage that applies"),
    "END": ("DIAG", False, False, "last stage (name or number)"),
    "EXISTING": ("skip", False, False, "skip | archive | delete: what to do with stages in the range that have results"),
    "SFM_SOURCE": ("opensfm", False, False, "opensfm (stages 1-4) | <finished ODM run dir> (read only; stages 1-4 not used)"),
    "OSFM_CONFIG": ("${REPO}/tools/opensfm_config.yaml", False, True, "OpenSfM config template"),
    "OSFM_PROCESSES": ("auto", False, False, "OpenSfM processes; auto = CPUs - 2"),
    "CAMERA_OVERRIDES": ("auto", False, True, "auto (brown from the EXIF focal) | OpenSfM camera_models_overrides.json"),
    "GSD": ("auto", False, False, "ortho cell (m); auto = flying height / focal"),
    "REF_CELL": ("2", False, False, "cell of the reference DSM from the SfM points (m)"),
    "REF_DSM": ("auto", False, True, "long-wave reference DSM: auto (the SfM DSM) | stereo (made in PAIRS) | <tif>"),
    "SHAPE_DSM": ("ref", False, True, "long-wave shape for BLEND: ref (= REF_DSM) | self (no long-wave change) | <tif>"),
    "BLEND_SIGMA": ("30", False, False, "BLEND long-wave scale (m)"),
    "BBA_ITERS": ("30", False, False, "BBA iterations"),
    "BBA_OPTIONS": ("", False, False, "extra bba.py options, e.g. --sig-gcp 0.03"),
    "TIE_WORKERS": ("2", False, False, "TIEPOINT worker processes"),
    "USE_TRACKS": ("no", False, False, "yes | no: run TRACKS and use them in BBA"),
    "PRODUCT_NAME": ("${BLOCK}_kjr", False, False, "product file prefix: <name>_dsm.tif, <name>_ortho.tif, ..."),
    "EXPORT_DIR": ("", False, True, "EXPORT: delivery root (files go to EXPORT_DIR/<BLOCK>)"),
    "REGIONAL_DEM": ("", False, True, "EXPORT: regional DEM for the context map of all blocks in EXPORT_DIR (empty: no map)"),
    "KJR_ENV": ("", False, True, "python env for the KJR stages (default: env.sh)"),
    "UAV_LIB": ("", False, True, "UAV matcher library (default: env.sh)"),
    "OSFM_ENV": ("/work/eridu_drone/kjr/opensfm/env", False, True, "conda env with OpenSfM built (tools/build_opensfm.sh)"),
    "OSFM_SRC": ("/work/eridu_drone/kjr/opensfm/src", False, True, "OpenSfM source checkout"),
}
NOISE = re.compile(r"Warn|warn|census variant|exact streaming|FutureWarning|diags_array|warnings.warn|^Note: In Python")


class CfgError(Exception):
    pass


# ---------------------------------------------------------------- control file
def read_cfg(path, seen=()):
    """-> {KEY: (raw value, dir of the file it came from)}"""
    path = os.path.abspath(path)
    if path in seen:
        raise CfgError(f"{path}: INCLUDE loop")
    if not os.path.isfile(path):
        raise CfgError(f"{path}: no such control file")
    here, out, inc_vals = os.path.dirname(path), {}, {}
    for n, line in enumerate(open(path, encoding="utf-8"), 1):
        line = re.sub(r"(^|\s)#.*", "", line).strip()
        if not line:
            continue
        if ":" not in line:
            raise CfgError(f"{path}:{n}: expected 'KEY : value': {line!r}")
        k, v = (x.strip() for x in line.split(":", 1))
        k = k.upper()
        if k == "INCLUDE":
            inc = v if os.path.isabs(v) else os.path.join(here, v)
            inc_vals.update(read_cfg(inc, seen + (path,)))
            continue
        if k not in KEYS:
            close = difflib.get_close_matches(k, list(KEYS) + ["INCLUDE"], 3, 0.6)
            raise CfgError(f"{path}:{n}: unknown key {k}" + (f" (did you mean {', '.join(close)}?)" if close else ""))
        if k in out:
            raise CfgError(f"{path}:{n}: {k} given twice")
        out[k] = (v, here)
    return {**inc_vals, **out}   # the including file wins, wherever the INCLUDE line is


def resolve(raw, overrides):
    """raw values + defaults -> {KEY: value}, ${KEY} substituted, paths absolute"""
    vals, base = {"REPO": REPO}, {}
    for k, (d, _, _, _) in KEYS.items():
        v, here = raw.get(k, ("", None))
        if k in overrides and overrides[k] is not None:
            v, here = str(overrides[k]), os.getcwd()
        vals[k], base[k] = (v, here) if v != "" else (d, None)
    for _ in range(10):
        changed = False
        for k, v in vals.items():
            nv = re.sub(r"\$\{(\w+)\}", lambda m: vals.get(m.group(1).upper(), m.group(0)), v)
            if nv != v:
                vals[k], changed = nv, True
        if not changed:
            break
    for k, v in vals.items():
        if "${" in v:
            raise CfgError(f"{k}: cannot resolve {v}")
    for k, (_, req, is_path, _) in KEYS.items():
        v = vals[k]
        if req and not v:
            raise CfgError(f"{k} is required")
        keyword = v.lower() in ("auto", "stereo", "ref", "self", "opensfm")
        if (is_path or k == "SFM_SOURCE") and v and not keyword and not os.path.isabs(v):
            vals[k] = os.path.normpath(os.path.join(base[k] or os.getcwd(), v))
    for k, allowed in (("EXISTING", ("skip", "archive", "delete")), ("USE_TRACKS", ("yes", "no"))):
        vals[k] = vals[k].lower()
        if vals[k] not in allowed:
            raise CfgError(f"{k} : {vals[k]} (allowed: {' | '.join(allowed)})")
    for k in ("REF_CELL", "BLEND_SIGMA", "BBA_ITERS", "TIE_WORKERS") + (("GSD",) if vals["GSD"] != "auto" else ()):
        try:
            float(vals[k])
        except ValueError:
            raise CfgError(f"{k} : {vals[k]} is not a number")
    del vals["REPO"]
    return vals


def stage_index(s):
    s = s.strip().upper()
    if s.isdigit() and 1 <= int(s) <= len(STAGES):
        return int(s) - 1
    if s in NAMES:
        return NAMES.index(s)
    raise CfgError(f"unknown stage {s!r} (use --stages)")


# ---------------------------------------------------------------- one block
class Job:
    def __init__(self, cfg_path, v, dry):
        self.cfg_path, self.v, self.dry = cfg_path, v, dry
        self.B, self.W, self.OUT = v["BLOCK"], v["WORK_DIR"], v["OUTPUT_DIR"]
        self.use_osfm = v["SFM_SOURCE"].lower() == "opensfm"
        self.ODM = os.path.join(self.W, "sfm/odm") if self.use_osfm else v["SFM_SOURCE"]
        self.D = os.path.join(self.W, "sfm/osfm")
        self.P = v["PRODUCT_NAME"]
        r = v["REF_DSM"]
        self.ref_stereo = r.lower() == "stereo"
        self.REF = (os.path.join(self.W, "ref/dsm_fused.tif") if self.ref_stereo else
                    os.path.join(self.ODM, "odm_dem/dsm.tif") if r.lower() == "auto" else r)
        s = v["SHAPE_DSM"]
        self.SHAPE = (self.REF if s.lower() == "ref" else
                      os.path.join(self.W, "stereo/dsm_fused.tif") if s.lower() == "self" else s)
        # centres the stereo search: the reference DSM, or the SfM DSM while the reference is made
        self.PRIOR = os.path.join(self.ODM, "odm_dem/dsm.tif") if self.ref_stereo else self.REF
        self.jobdir = os.path.join(self.W, "job")
        self.stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # which stages apply, and their outputs (paths that a rerun archives / deletes)
    def applies(self, s):
        if s.startswith("SFM_") or s == "GEOREF":
            return self.use_osfm
        if s == "TRACKS":
            return self.v["USE_TRACKS"] == "yes"
        if s == "EXPORT":
            return bool(self.v["EXPORT_DIR"])
        return True

    def outputs(self, s):
        D, W, O, P = self.D, self.W, self.OUT, self.P
        d = {
            "SFM_FEATURE": [f"{D}/exif", f"{D}/camera_models.json", f"{D}/camera_models_overrides.json", f"{D}/features",
                            f"{D}/reports/features", f"{D}/reports/features.json"],
            "SFM_MATCH": [f"{D}/reference_lla.json", f"{D}/matches", f"{D}/tracks.csv", f"{D}/reports/matches.json",
                          f"{D}/reports/tracks.json"],
            "SFM_RECON": [f"{D}/reconstruction.json", f"{D}/reports/reconstruction.json"],
            "GEOREF": [self.ODM] if self.use_osfm else [],
            "PAIRS": [f"{W}/ref"],
            "TIEPOINT": [f"{W}/ties", f"{W}/sift"],
            "TRACKS": [f"{W}/tracks.npz"],
            "BBA": [f"{W}/bba"],
            "UNDISTORT": [f"{W}/images"],
            "DSM_EXTRACTION": [f"{W}/stereo/pairs", f"{W}/stereo/stereo_meta.json", f"{W}/stereo/pairs.log"],
            "FUSE": [f"{W}/stereo/{f}" for f in ("dsm_fused.tif", "dsm_raw.tif", "npairs.tif", "pair_planes.json")],
            "BLEND": [f"{W}/final/dsm_fused.tif"],
            "ORTHO": [f"{W}/final/ortho.tif"],
            "PRODUCTS": [f"{O}/{P}_{x}" for x in ("ortho.tif", "dsm.tif", "ortho_envi.tif", "dsm_envi.tif", "ortho.jpg",
                         "ortho.jpg.aux.xml", "ortho.wld", "ortho_85pct.jpg", "ortho_85pct.jpg.aux.xml",
                         "ortho_85pct.wld", "dsm_hillshade.jpg", "dsm_hillshade.jpg.aux.xml", "dsm_hillshade.wld")],
            "DIAG": [f"{W}/diag", f"{O}/{P}_dsm_color_hillshade_small.jpg"],
            "EXPORT": [os.path.join(self.v["EXPORT_DIR"], self.B)] if self.v["EXPORT_DIR"] else [],
        }
        return d[s]

    def has_results(self, s):
        """outputs present (for runs made before job control: no marker)"""
        W, D = self.W, self.D
        key = {
            "SFM_FEATURE": [f"{D}/features"], "SFM_MATCH": [f"{D}/tracks.csv"],
            "SFM_RECON": [f"{D}/reconstruction.json"], "GEOREF": [f"{self.ODM}/odm_dem/dsm.tif"],
            "PAIRS": [f"{W}/ref/stereo_meta.json"] + ([f"{W}/ref/dsm_fused.tif"] if self.ref_stereo else []),
            "TRACKS": [f"{W}/tracks.npz"], "BBA": [f"{W}/bba/bba.json"],
            "FUSE": [f"{W}/stereo/dsm_fused.tif"], "BLEND": [f"{W}/final/dsm_fused.tif"],
            "ORTHO": [f"{W}/final/ortho.tif"], "PRODUCTS": [f"{self.OUT}/{self.P}_dsm.tif", f"{self.OUT}/{self.P}_ortho.jpg"],
            "DIAG": [f"{W}/diag/imgfixed.npy"],
        }
        if s in key:
            return all(os.path.exists(p) for p in key[s])
        if s == "TIEPOINT":
            m = f"{W}/ref/stereo_meta.json"
            return os.path.isfile(m) and len([f for f in os.listdir(f"{W}/ties") if f.endswith(".npz")]
                                             if os.path.isdir(f"{W}/ties") else []) >= len(json.load(open(m))["pairs"])
        if s == "UNDISTORT":
            n = len(self.images())
            return os.path.isdir(f"{W}/images") and len([f for f in os.listdir(f"{W}/images") if f.endswith(".tif")]) >= n
        if s == "DSM_EXTRACTION":   # pairs that fail write nothing: a count cannot tell; FUSE ran = complete
            return os.path.exists(f"{W}/stereo/dsm_fused.tif")
        return False                # EXPORT: marker only

    def images(self):
        return sorted(f for f in os.listdir(self.v["INPUT_DIR"]) if f.lower().endswith((".jpg", ".jpeg")))

    # markers
    def marker(self, s):
        return os.path.join(self.jobdir, "done", f"{NAMES.index(s) + 1:02d}_{s}")

    def stale(self):
        f = os.path.join(self.jobdir, "stale")
        return set(open(f).read().split()) if os.path.exists(f) else set()

    def set_stale(self, st):
        os.makedirs(self.jobdir, exist_ok=True)
        open(os.path.join(self.jobdir, "stale"), "w").write("\n".join(sorted(st, key=NAMES.index)) + ("\n" if st else ""))

    def complete(self, s):
        return s not in self.stale() and (os.path.exists(self.marker(s)) or self.has_results(s))

    # logging and running
    def say(self, msg):
        line = f"== {datetime.datetime.now():%F %T} {self.B} {msg}"
        print(line, flush=True)
        if not self.dry:
            open(os.path.join(self.W, "kjr_run.log"), "a").write(line + "\n")

    def run(self, s, cmd, env=None, cwd=None, filt=True):
        """run one command of stage s; full output in job/logs/<nn>_<stage>.log, stop on failure"""
        shown = " ".join(str(c) for c in cmd)
        if self.dry:
            print(f"     $ {shown if len(shown) < 400 else shown[:400] + ' ...'}")
            return
        log = os.path.join(self.jobdir, "logs", f"{NAMES.index(s) + 1:02d}_{s}.log")
        with open(log, "a") as lf:
            lf.write(f"\n$ {shown}\n"); lf.flush()
            p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                 env=env or self.env, cwd=cwd, bufsize=1)
            tail = []
            for line in p.stdout:
                if filt and NOISE.search(line):
                    continue
                lf.write(line)
                tail = (tail + [line.rstrip()])[-25:]
            rc = p.wait()
        if rc != 0:
            self.say(f"{s} FAILED ({rc}): {shown}\n" + "\n".join("   | " + t for t in tail) + f"\n   full log: {log}")
            raise SystemExit(f"{self.B}: {s} failed")
        return tail

    def py(self, script):
        return [self.PY, os.path.join(self.code, script)]

    # ------------------------------------------------------------ setup
    def check_inputs(self):
        v = self.v
        for k in ("INPUT_DIR",):
            if not os.path.isdir(v[k]):
                raise CfgError(f"{k}: no such directory {v[k]}")
        if not self.images():
            raise CfgError(f"INPUT_DIR: no *.jpg in {v['INPUT_DIR']}")
        if not os.path.isfile(v["GCP_FILE"]):
            raise CfgError(f"GCP_FILE: no such file {v['GCP_FILE']}")
        if not self.use_osfm:
            for f in ("opensfm/reconstruction.json", "opensfm/undistorted/reconstruction.json",
                      "odm_georeferencing/coords.txt", "odm_orthophoto/odm_orthophoto.tif"):
                if not os.path.exists(os.path.join(self.ODM, f)):
                    raise CfgError(f"SFM_SOURCE: {self.ODM} is not a finished ODM run (no {f})")
        for k in ("REF_DSM", "SHAPE_DSM", "REGIONAL_DEM"):
            if v[k] and v[k].lower() not in ("auto", "stereo", "ref", "self") and not os.path.isfile(v[k]):
                raise CfgError(f"{k}: no such file {v[k]}")
        if v["CAMERA_OVERRIDES"] != "auto" and not os.path.isfile(v["CAMERA_OVERRIDES"]):
            raise CfgError(f"CAMERA_OVERRIDES: no such file {v['CAMERA_OVERRIDES']}")
        if self.use_osfm and not os.path.isfile(os.path.join(v["OSFM_SRC"], "bin/opensfm_main.py")):
            raise CfgError(f"OSFM_SRC: no bin/opensfm_main.py in {v['OSFM_SRC']} (tools/build_opensfm.sh)")

    def snapshot(self):
        """copy the repo (code only) into the job dir and run everything from the copy"""
        jd = os.path.join(self.jobdir, self.stamp)
        self.code = os.path.join(jd, "code")
        if self.dry:
            self.code = REPO
        else:
            os.makedirs(os.path.join(self.jobdir, "logs"), exist_ok=True)
            os.makedirs(os.path.join(self.jobdir, "done"), exist_ok=True)
            shutil.copytree(REPO, self.code, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc"))
            shutil.copy(self.cfg_path, os.path.join(jd, "control.cfg"))
            with open(os.path.join(jd, "resolved.cfg"), "w") as f:
                f.write(f"# resolved from {self.cfg_path} at {self.stamp}\n")
                for k, x in self.v.items():
                    f.write(f"{k:<16} : {x}\n")
            g = lambda *a: subprocess.run(["git", "-C", REPO, *a], capture_output=True, text=True).stdout
            open(os.path.join(jd, "git.txt"), "w").write(g("log", "-1", "--format=%H %s") + g("status", "--short")
                                                         + g("diff"))
        # environment of env.sh (KJR_ENV, UAV_LIB from the control file if given), from the copy
        e = dict(os.environ)
        for k in ("KJR_ENV", "UAV_LIB"):
            if self.v[k]:
                e[k] = self.v[k]
        out = subprocess.run(["bash", "-c", f'source "{self.code}/env.sh" && env -0'], env=e, capture_output=True,
                             check=True).stdout
        self.env = dict(x.split("=", 1) for x in out.decode().split("\0") if "=" in x)
        self.PY = self.env["PY"]

    # ------------------------------------------------------------ stages
    def osfm_setup(self):
        D, v = self.D, self.v
        if self.dry:
            return
        os.makedirs(f"{D}/images", exist_ok=True)
        for f in self.images() + [f for f in os.listdir(v["INPUT_DIR"]) if f.lower().endswith((".tif", ".tiff", ".png"))]:
            dst = f"{D}/images/{f}"
            if not os.path.lexists(dst):
                os.symlink(os.path.join(v["INPUT_DIR"], f), dst)
        shutil.copy(v["GCP_FILE"], f"{D}/gcp_list.txt")
        n = os.cpu_count() or 3
        proc = max(n - 2, 1) if v["OSFM_PROCESSES"] == "auto" else int(v["OSFM_PROCESSES"])
        with open(f"{D}/config.yaml", "w") as f:
            f.write(open(v["OSFM_CONFIG"]).read().rstrip("\n") + f"\nprocesses: {proc}\n")

    def osfm(self, s, step):
        e = dict(os.environ)   # not the KJR env: OpenSfM has its own python and libraries
        e["PATH"] = os.path.join(self.v["OSFM_ENV"], "bin") + ":" + e.get("PATH", "")
        e.pop("PYTHONPATH", None)
        self.run(s, [os.path.join(self.v["OSFM_ENV"], "bin/python"), os.path.join(self.v["OSFM_SRC"], "bin/opensfm_main.py"),
                     step, self.D], env=e, filt=False)

    def stage(self, s):
        v, W, ODM, D = self.v, self.W, self.ODM, self.D
        if s == "SFM_FEATURE":
            self.osfm_setup()
            if v["CAMERA_OVERRIDES"] != "auto":
                if not self.dry:
                    shutil.copy(v["CAMERA_OVERRIDES"], f"{D}/camera_models_overrides.json")
            else:   # every EXIF camera as brown (k1-k3, p1, p2) with its EXIF focal, self-calibrated
                self.osfm(s, "extract_metadata")
                if not self.dry:
                    cams = json.load(open(f"{D}/camera_models.json"))
                    o = {k: {"projection_type": "brown", "width": c["width"], "height": c["height"],
                             "focal_x": c.get("focal_x", c.get("focal")), "focal_y": c.get("focal_x", c.get("focal")),
                             "c_x": 0.0, "c_y": 0.0, "k1": 0.0, "k2": 0.0, "k3": 0.0, "p1": 0.0, "p2": 0.0}
                         for k, c in cams.items()}
                    json.dump(o, open(f"{D}/camera_models_overrides.json", "w"), indent=1)
                    self.say("cameras: " + ", ".join(f"{k} focal {c['focal_x']:.4f}" for k, c in o.items()))
            self.osfm(s, "extract_metadata")
            self.osfm(s, "detect_features")
        elif s == "SFM_MATCH":
            self.osfm_setup()
            self.osfm(s, "match_features")
            self.osfm(s, "create_tracks")
        elif s == "SFM_RECON":
            self.osfm_setup()
            self.osfm(s, "reconstruct")
        elif s == "GEOREF":
            t = self.run(s, self.py("kjr/osfm_import.py") + [D, ODM, "--cell", v["REF_CELL"]]
                         + ([] if v["GSD"] == "auto" else ["--gsd", v["GSD"]]))
            for line in t or []:
                if re.match(r"reconstruction:|WARNING|ENU|\d+ shots, f ", line):
                    self.say("   " + line)
        elif s == "PAIRS":
            if self.ref_stereo:
                self.run(s, self.py("kjr/stereo.py") + [ODM, f"{W}/ref", "--prior-dsm", self.PRIOR])
                self.run(s, self.py("kjr/fuse.py") + [f"{W}/ref"])
            else:
                self.run(s, self.py("kjr/stereo.py") + [ODM, f"{W}/ref", "--meta-only", "--prior-dsm", self.PRIOR])
        elif s == "TIEPOINT":
            self.run(s, self.py("kjr/ties.py") + [f"{W}/ref/stereo_meta.json", v["INPUT_DIR"], W,
                                                  "--workers", v["TIE_WORKERS"], "--keep-sift"])
            if v["USE_TRACKS"] != "yes" and not self.dry:
                shutil.rmtree(f"{W}/sift", ignore_errors=True)   # SIFT cache (~80 MB per image), TRACKS reuses it
        elif s == "TRACKS":
            self.run(s, self.py("kjr/tracks.py") + [f"{W}/ref/stereo_meta.json", v["INPUT_DIR"], W])
            if not self.dry:
                shutil.rmtree(f"{W}/sift", ignore_errors=True)
        elif s == "BBA":
            trk = ["--tracks", f"{W}/tracks.npz"] if v["USE_TRACKS"] == "yes" else []
            self.run(s, self.py("kjr/bba.py") + [ODM, f"{W}/ties", v["GCP_FILE"], f"{W}/bba", "--iters", v["BBA_ITERS"]]
                     + trk + v["BBA_OPTIONS"].split())
            if not self.dry:
                lines = open(f"{W}/bba/bba.log").read().splitlines()
                its = [x for x in lines if x.startswith("it ")]
                for line in its[-1:] + [x for x in lines if re.search(r"GCP . triangulated|grid correction", x)][-6:]:
                    self.say("   " + line)
        elif s == "UNDISTORT":
            self.run(s, self.py("kjr/undistort.py") + [f"{W}/bba/odm", v["INPUT_DIR"], f"{W}/images"]
                     + [f + ".tif" for f in self.images()])
        elif s == "DSM_EXTRACTION":
            self.run(s, self.py("kjr/stereo.py") + [f"{W}/bba/odm", f"{W}/stereo", "--img-dir", f"{W}/images",
                                                    "--prior-dsm", self.PRIOR])
        elif s == "FUSE":
            self.run(s, self.py("kjr/fuse.py") + [f"{W}/stereo"])
        elif s == "BLEND":
            if not self.dry:
                os.makedirs(f"{W}/final", exist_ok=True)
            self.say(f"   reference {self.SHAPE}, sigma {v['BLEND_SIGMA']} m")
            self.run(s, self.py("kjr/lowblend.py") + [f"{W}/stereo/dsm_fused.tif", self.SHAPE, f"{W}/final/dsm_fused.tif",
                                                      "--sigma", v["BLEND_SIGMA"], "--gcp", v["GCP_FILE"]])
        elif s == "ORTHO":
            self.run(s, self.py("kjr/ortho.py") + [f"{W}/bba/odm", f"{W}/final", "--img-dir", f"{W}/images"])
        elif s == "PRODUCTS":
            self.run(s, [os.path.join(self.code, "kjr/make_products.sh"), self.B, f"{W}/final", self.OUT, self.P])
        elif s == "DIAG":
            if not self.dry:
                os.makedirs(f"{W}/diag", exist_ok=True)
            self.run(s, self.py("tools/color_hillshade.py") + [f"{W}/final/dsm_fused.tif",
                                                               f"{self.OUT}/{self.P}_dsm_color_hillshade_small.jpg"])
            t = self.run(s, self.py("tools/gcp_check.py") + [v["GCP_FILE"], self.SHAPE, f"{W}/stereo/dsm_fused.tif",
                                                             f"{W}/final/dsm_fused.tif"])
            for line in t or []:
                self.say("   " + line)
            t = self.run(s, self.py("tools/diag_imgfixed.py") + [f"{W}/final/dsm_fused.tif", f"{W}/bba/odm",
                                                                 f"{W}/diag/imgfixed.npy"])
            for line in t or []:
                if "hp DSM" in line:
                    self.say("   " + line)
        elif s == "EXPORT":
            X, E = os.path.join(v["EXPORT_DIR"], self.B), v["EXPORT_DIR"]
            if not self.dry:
                os.makedirs(X, exist_ok=True)
                for f in (f"{self.P}_dsm.tif", f"{self.P}_ortho.tif", f"{self.P}_dsm_color_hillshade_small.jpg"):
                    shutil.copy(os.path.join(self.OUT, f), X)
            self.run(s, self.py("tools/small_ortho.py") + [f"{self.OUT}/{self.P}_ortho.tif",
                                                           f"{self.OUT}/{self.P}_dsm_color_hillshade_small.jpg",
                                                           f"{X}/{self.P}_ortho_small.jpg"])
            self.run(s, self.py("tools/image_footprints.py") + [self.P, W, self.OUT, f"{X}/{self.B}_image_footprints.png",
                                                                f"{X}/{self.P}_ortho_small.jpg"])
            if v["REGIONAL_DEM"]:   # context map of every block exported so far
                blocks = []
                for b in sorted(os.listdir(E)) if os.path.isdir(E) else []:
                    dsm = [f for f in os.listdir(os.path.join(E, b)) if f.endswith("_dsm.tif")] \
                        if os.path.isdir(os.path.join(E, b)) else []
                    if dsm:
                        blocks.append(f"{b}={os.path.join(E, b, dsm[0])}")
                if self.dry:
                    blocks = blocks or [f"{self.B}=<{X}/{self.P}_dsm.tif>"]
                self.run(s, self.py("tools/context_map.py") + [v["REGIONAL_DEM"], f"{E}/uav_context_map.png"] + blocks)

    # ------------------------------------------------------------ driver
    def clear(self, s, how):
        """move (archive) or delete the old outputs of stage s"""
        old = [p for p in self.outputs(s) if os.path.lexists(p)]
        if not self.dry and os.path.exists(self.marker(s)):
            os.remove(self.marker(s))
        if not old:
            return
        def rel(p):   # path inside the archive
            for root, sub in ((self.W, ""), (self.OUT, "products"), (self.v["EXPORT_DIR"], "export")):
                if root and p.startswith(root.rstrip("/") + "/"):
                    return os.path.join(sub, os.path.relpath(p, root))
            return os.path.basename(p)
        if how == "archive":
            dst = os.path.join(self.W, "old", self.stamp)
            self.say(f"{s}: old results -> {dst}")
        else:
            self.say(f"{s}: deleting old results: " + ", ".join(rel(p) for p in old))
        if self.dry:
            return
        for p in old:
            if how == "archive":
                d = os.path.join(dst, rel(p))
                os.makedirs(os.path.dirname(d), exist_ok=True)
                shutil.move(p, d)
            elif os.path.isdir(p) and not os.path.islink(p):
                shutil.rmtree(p)
            else:
                os.remove(p)

    def go(self, start, end):
        v = self.v
        a = stage_index(start) if start else next(i for i, s in enumerate(NAMES) if self.applies(s))
        b = stage_index(end or v["END"])
        if b < a:
            raise CfgError(f"END {NAMES[b]} comes before START {NAMES[a]}")
        if not self.applies(NAMES[a]):
            why = {"TRACKS": "USE_TRACKS : no", "EXPORT": "EXPORT_DIR is empty"}.get(NAMES[a], f"SFM_SOURCE is {v['SFM_SOURCE']}")
            raise CfgError(f"START {NAMES[a]} does not apply ({why})")
        if NAMES[b] == "EXPORT" and not self.applies("EXPORT"):
            raise CfgError("END EXPORT needs EXPORT_DIR")
        self.check_inputs()
        rng = [s for s in NAMES[a:b + 1] if self.applies(s)]
        unused = [s for s in NAMES[a:b + 1] if not self.applies(s)]
        # stages before START must be complete and not stale
        st = self.stale()
        for s in NAMES[:a]:
            if not self.applies(s):
                continue
            if s in st:
                raise CfgError(f"{s} (before START) is stale: an earlier stage was run again after it. "
                               f"Set START to {s} or earlier.")
            if not self.complete(s):
                raise CfgError(f"{s} (before START) has no results in {self.W}. Set START to {s} or earlier.")
        self.snapshot()
        if not self.dry:
            os.makedirs(self.OUT, exist_ok=True)
        self.say(f"start {v['TITLE'] or ''} | {self.cfg_path} | stages {rng[0]} .. {rng[-1]} | EXISTING {v['EXISTING']}"
                 f" | code {os.path.relpath(self.code, self.W) if not self.dry else REPO}"
                 + (f" | not used: {' '.join(unused)}" if unused else ""))
        for s in rng:
            n = NAMES.index(s) + 1
            if v["EXISTING"] != "skip":
                self.clear(s, v["EXISTING"])
            elif s in self.stale():   # results made from older inputs: invalid, and the tools resume from them
                if any(os.path.lexists(p) for p in self.outputs(s)):
                    self.say(f"{n:2d} {s}: results are stale (an earlier stage ran again)")
                self.clear(s, "delete")
            elif self.complete(s):
                self.say(f"{n:2d} {s}: complete, skipped")
                continue
            self.say(f"{n:2d} {s}: {dict(STAGES)[s]}")
            t0 = datetime.datetime.now()
            self.stage(s)
            if not self.dry:
                st = (self.stale() - {s}) | {x for x in NAMES[n:] if self.applies(x)}
                self.set_stale(st)
                open(self.marker(s), "w").write(f"{datetime.datetime.now():%F %T} {self.stamp}\n")
            self.say(f"{n:2d} {s}: done ({str(datetime.datetime.now() - t0).split('.')[0]})")
        later = sorted((s for s in self.stale() if NAMES.index(s) > b), key=NAMES.index)
        if later and not self.dry:
            self.say("note: results of later stages are now stale (run them again): " + " ".join(later))
        self.say("finished")


def main():
    ap = argparse.ArgumentParser(description="uav_kjr: run blocks from control files (*.cfg)",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n", 2)[2])
    ap.add_argument("cfg", nargs="*", help="control files; several = one block after the other")
    ap.add_argument("--start", help="override START")
    ap.add_argument("--end", help="override END")
    ap.add_argument("--existing", choices=("skip", "archive", "delete"), help="override EXISTING")
    ap.add_argument("--dry-run", action="store_true", help="check, list the stages and commands, run nothing")
    ap.add_argument("--stages", action="store_true", help="list the stages")
    ap.add_argument("--keys", action="store_true", help="list the control file keys")
    a = ap.parse_args()
    if a.stages:
        for i, (s, d) in enumerate(STAGES, 1):
            print(f"{i:2d} {s:<15} {d}")
        return
    if a.keys:
        for k, (d, req, _, h) in KEYS.items():
            print(f"{k:<17} {'(required)' if req else d or '(empty)':<36} {h}")
        return
    if not a.cfg:
        ap.error("no control file")
    failed = []
    for c in a.cfg:
        try:
            v = resolve(read_cfg(c), {"EXISTING": a.existing})
            Job(os.path.abspath(c), v, a.dry_run).go(a.start or v["START"], a.end)
        except CfgError as e:
            print(f"{c}: {e}", file=sys.stderr)
            failed.append(c)
        except SystemExit as e:
            print(e, file=sys.stderr)
            failed.append(c)
    if len(a.cfg) > 1:
        print(f"== {len(a.cfg) - len(failed)} of {len(a.cfg)} blocks finished" + (f"; failed: {' '.join(failed)}" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
