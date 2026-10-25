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
