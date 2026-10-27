# uav_kjr (ODM-KJR)

DSM and orthophoto of a drone block from the raw images and a GCP list, without the in-DSM
artefacts of ODM: concentric ripples (camera-model error fixed to image coordinates) and steps
between image strips (pose error). **ODM is not needed** (licences: THIRD_PARTY_LICENSES.md).

upstream OpenSfM (BSD) SfM with the GCPs → sub-pixel ties (SIFT, MAGSAC++, Gruen ALSC) → own
bundle adjustment with radial + 2-D grid image correction → Lanczos undistortion → census-SGM
per pair → robust fusion → long-wave blend with the SfM reference DSM → orthophoto →
COG / ENVI / JPG products.

```bash
tools/build_opensfm.sh /path/opensfm          # once (upstream OpenSfM, unmodified)
./kjr_run cfg/example.cfg                     # one block from a control file
./kjr_run cfg/eridu/*.cfg                     # several blocks, one after the other
./kjr_run JJ3.cfg --start BBA --end DIAG --existing archive   # stage range / old results
./kjr_run --stages ; ./kjr_run --keys ; ./kjr_run JJ3.cfg --dry-run
```

Control file (`*.cfg`, "KEY : value"): `BLOCK`, `INPUT_DIR`, `GCP_FILE`, `WORK_DIR`, `OUTPUT_DIR`,
`START`/`END` (16 stages: SFM_FEATURE SFM_MATCH SFM_RECON GEOREF PAIRS TIEPOINT TRACKS BBA
UNDISTORT DSM_EXTRACTION FUSE BLEND ORTHO PRODUCTS DIAG EXPORT), `EXISTING : skip | archive |
delete`, options; `INCLUDE` for shared settings. Each run records the control file, all
resolved values, the git state and a copy of the code it ran in `<WORK_DIR>/job/`.
`SFM_SOURCE : <finished ODM run>` starts from ODM results instead of OpenSfM.
`run.sh` / `run_block.sh` remain as wrappers for the old command lines.

JJ-3 (117 images, 5 GCPs): OpenSfM 11 min, the BBA lands where it did from ODM (GCP residuals
within 4 cm), final DSM GCP height RMS 5.3 cm (9.9 cm from ODM).

User manual (Korean) with the standard templates: [manual/](manual/README.md)
(`manual/templates/new_block.sh` starts a new block). Technical manual (Korean): [MANUAL.md](MANUAL.md)

| | ODM cameras | ODM-KJR |
|---|---|---|
| image-locked ripple (std / even-odd corr) | 1.14 cm / 0.84 | 0.32 cm / 0.19 |
| steps between consecutive pairs, east strip | up to 18 cm | up to 8 cm |
| GCP height RMS | 3.8 cm | 5.2 cm |

KJ-1 (DJI FC6310, 68 images, 5 GCPs).

Speed (KJ-3, 61 images, 215 pairs, RTX 4060 8 GB): 1 h 47 min → 32 min per block, same result
(GCP RMS 7.1 cm, ripple 0.51 cm). GPU guided SIFT matching, a C++/OpenMP point-elimination
Schur kernel for the BBA (`kjr/csrc/schur.cpp`, built on first use), GPU fusion statistics,
threaded stereo and ortho. Details in MANUAL.md §2.1.

Deliverables: `tools/export_products.sh` (DSM, ortho, small hillshade / ortho, image
footprint map per block, context map on a regional DEM).

## Own SfM (`run_noodm.sh`, failed — do not use)

EXIF/XMP initial cameras, GCP similarity alignment, self-calibrating bundle adjustment and a
tie-point reference DSM replace the ODM run. It processes a block end to end, but with nadir-only
parallel strips and few GCPs the lens distortion and the block curvature are not separable:
KJ-3 leave-one-out check points RMS 80 cm self-calibrated, 37 cm with the camera calibration of
another block (`CAMERA=`), 21 cm starting from ODM. Give a camera calibration and check with
`CHECK=1`; see MANUAL 4.2. On JJ-3 it failed outright (focal 10 % off, metre-level
errors at the block ends); use `kjr_run` (OpenSfM) instead.
