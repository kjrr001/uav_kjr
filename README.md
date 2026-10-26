# uav_kjr (ODM-KJR)

Rebuilds the DSM and orthophoto of an ODM drone block without the in-DSM artefacts:
concentric ripples (camera-model error fixed to image coordinates) and steps between
image strips (pose error).

ODM solution → sub-pixel ties (SIFT, MAGSAC++, Gruen ALSC) → own bundle adjustment with
radial + 2-D grid image correction → Lanczos undistortion → census-SGM per pair → robust
fusion → long-wave blend with the ODM-camera DSM → orthophoto → COG / ENVI / JPG products.

```bash
./run_block.sh KJ1 <odm run> <image dir> <gcp_list.txt> <work dir> <products dir>
```

Manual (Korean): [MANUAL.md](MANUAL.md)

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

## Without ODM (`run_noodm.sh`, experimental)

EXIF/XMP initial cameras, GCP similarity alignment, self-calibrating bundle adjustment and a
tie-point reference DSM replace the ODM run. It processes a block end to end, but with nadir-only
parallel strips and few GCPs the lens distortion and the block curvature are not separable:
KJ-3 leave-one-out check points RMS 80 cm self-calibrated, 37 cm with the camera calibration of
another block (`CAMERA=`), 21 cm starting from ODM. Give a camera calibration and check with
`CHECK=1`; see MANUAL 4.1.
