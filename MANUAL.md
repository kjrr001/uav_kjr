# ODM-KJR 매뉴얼

ODM으로 처리한 드론 블록에서 DSM과 정사영상을 다시 만드는 도구다. ODM 결과의 카메라와 GCP 정합은 그대로 출발점으로 쓰고, 다음 세 단계를 자체 코드로 바꾼다.

1. **tie point**: 원본 JPEG에서 SIFT → MAGSAC++ → Gruen ALSC로 sub-pixel tie point를 뽑는다.
2. **번들 조정(BBA)**: 자체 번들 조정에 사진 좌표 보정(반경 보정 + 2차원 격자 보정)을 더한다.
3. **DSM**: 사진 쌍마다 census-SGM 정합을 한 뒤 강건하게 융합한다.

목표는 측지 정확도가 아니라 **보기 좋은 DSM**이다. 한 DSM 안의 동심원 무늬, 블록 사이 단차, 줄무늬를 없애는 것이 목적이다.
현장에 식생이 없다는 전제(DSM = DTM)로 만들었으므로 DTM은 만들지 않는다.

---

## 1. 왜 이렇게 만들었나 (KJ-1에서 확인한 원인)

| 증상 | 원인 | 대응 |
|---|---|---|
| DSM의 동심원·물결 무늬 (1–3 cm) | **카메라 모델 오차**. 무늬가 지면이 아니라 사진 좌표에 고정되어 있다. 모든 사진을 사진 좌표로 겹쳐 평균하면 짝수 사진 묶음과 홀수 사진 묶음의 무늬가 상관 0.8–0.9로 일치한다. ODM 자체 DSM(OpenMVS)에도 같은 무늬가 있으므로 정합기 문제가 아니다. | 번들 조정에서 사진 좌표 보정을 추정한다: 반경 30 px 간격 + 2차원 격자 64 px 간격 |
| 동쪽 줄의 연속 사진 쌍 사이 단차 (10–18 cm) | ODM 카메라 자세 오차 | 자체 번들 조정 (tie point 약 108만 점) |
| 고리 모양의 정합 편향 | OpenSfM 왜곡 보정의 bilinear 보간 | Lanczos로 왜곡 보정한 사진을 따로 만든다 |
| 번들 조정 후 블록 전체가 휨 (±30 cm) | GCP 5점에 비해 보정 자유도가 크다 | 30 m보다 긴 파장은 ODM 카메라로 만든 참조 DSM에서 가져온다 (`lowblend`) |

시도했다가 버린 방법:
- **내부표정(f, k1–k3)을 풀어서 추정**: 1 m 크기의 돔이 생겼다. 그래서 내부표정은 ODM 값으로 고정한다.
- **사진 쌍별 평면 보정, 쌍 전체 동시 보정**: 오히려 단차를 만들었다. 그래서 `fuse.py`의 기본값을 `--iter 0`으로 둔다.

---

## 2. 처리 단계

`run_block.sh` 하나가 아래 단계를 순서대로 실행한다. 단계마다 결과가 이미 있으면 건너뛴다. 그래서 중간에 실패해도 다시 실행하면 이어서 진행한다.

| 단계 | 하는 일 | 결과 (`<work>/`) | KJ-1 소요 (RTX 4060 8 GB) |
|---|---|---|---|
| `ref` | ODM 카메라로 stereo + 융합 → 사진 쌍 목록, 긴 파장 참조 DSM | `ref/stereo_meta.json`, `ref/dsm_fused.tif` | 약 45분 |
| `ties` | 원본 JPEG에서 tie point 추출 (GPU ALSC) | `ties/<a>__<b>.npz` | 약 50분 (242쌍) |
| `bba` | 번들 조정 (LM + Schur 소거, Cauchy 가중) | `bba/bba.json`, `bba/bba.log`, `bba/odm/` | 약 25분 (30회 반복) |
| `undist` | 새 카메라로 Lanczos 왜곡 보정 | `images/*.JPG.tif` | 약 1분 |
| `stereo` | 사진 쌍별 census-SGM (pfsgbm, MGM, CUDA) | `stereo/pairs/*.npz` | 약 25분 |
| `fuse` | 강건 융합 (중앙값 + NMAD, 가장자리·외곽 가중) | `stereo/dsm_fused.tif` | 약 20분 |
| `blend` | 긴 파장(σ 30 m)은 `ref`, 세부는 `stereo` | `final/dsm_fused.tif` | 1분 |
| `ortho` | 최종 DSM으로 정사영상 (가장 수직에 가까운 사진에 가중) | `final/ortho.tif` | 약 5분 |
| `products` | COG, ENVI, JPG 제품 파일 | `<products>/` | 약 3분 |
| `diag` | 컬러 힐쉐이드, GCP 높이, 사진 좌표 무늬 검사 | `diag/`, 로그 | 약 2분 |

블록 하나에 모두 **약 3시간** 걸린다. 디스크는 작업 폴더에 약 10 GB, 제품 파일에 약 11 GB가 든다.

---

## 3. 설치

필요한 것:
- Linux, NVIDIA GPU (8 GB면 충분), CUDA 드라이버
- Python 환경 (conda prefix): numpy, scipy, opencv, gdal, matplotlib, pillow, torch, cuda-python
- UAV 정합 라이브러리 `UAV_LIB`: `tiepoints`(SIFT + MAGSAC++ + ALSC)와 `matcher`(census SGM). 현재 위치는 `/work10/planet_stereo/UAV`이다.
- GDAL 명령행 도구: `gdal_translate`, `gdalwarp`, `gdaldem`
- 끝난 ODM GCP 처리 폴더. 다음 파일이 있어야 한다:
  - `opensfm/reconstruction.json`
  - `opensfm/undistorted/reconstruction.json` 및 `opensfm/undistorted/images/`
  - `odm_georeferencing/coords.txt`
  - `odm_orthophoto/odm_orthophoto.tif`

`env.sh`가 경로를 정한다. 경로를 바꾸려면 실행 전에 환경변수로 지정한다.

```bash
export KJR_ENV=/path/to/conda/env        # 기본 /work/eridu_drone/kjr/env
export UAV_LIB=/path/to/UAV              # 기본 /work10/planet_stereo/UAV
```

> 드라이버가 CUDA 12.2(535)라면 torch cu121 휠의 NVRTC 12.1을 써야 한다. 더 새로운 NVRTC가 만든 PTX는 드라이버가 거부한다. `env.sh`가 `LD_LIBRARY_PATH`를 그렇게 맞춘다.

---

## 4. 사용법

```bash
./run_block.sh <블록> <ODM 처리 폴더> <원본 사진 폴더> <GCP 목록> <작업 폴더> <제품 폴더>
```

KJ-1 예:

```bash
cd /work/eridu_drone/kjr/odm-kjr
./run_block.sh KJ1 \
    /work/eridu_drone/odm_data/KJ1_gcp \
    /work/eridu_drone/KJ-1 \
    /work/eridu_drone/final_output_kj1/gcp_list_final.txt \
    /work/eridu_drone/kjr/runs/KJ1 \
    /work/eridu_drone/kjr/products/KJ1
```

오래 걸리므로 `nohup ... &`로 띄우고 `<작업 폴더>/run_block.log`를 본다.

**ODM 폴더와 원본 사진은 읽기만 한다.** 모든 결과는 작업 폴더와 제품 폴더에만 쓴다.

### 환경변수 옵션

| 변수 | 기본값 | 의미 |
|---|---|---|
| `STAGES` | 전체 | 실행할 단계만 지정. 예: `STAGES="blend ortho products"` |
| `REF_DSM` | (없음) | 긴 파장 참조 DSM을 직접 지정한다. 그러면 `ref` 단계는 사진 쌍 목록만 만들고 약 45분을 아낀다. |
| `SIGMA` | 30 | `blend`의 Gaussian σ (m) |
| `BBA_ITERS` | 30 | 번들 조정 반복 횟수 |
| `TIE_WORKERS` | 3 | tie point 병렬 작업 수. GPU 메모리가 모자라면 자동으로 재시도한다. |

### 다시 돌리기

- 단계마다 결과 파일이 있으면 건너뛴다. 그 단계를 다시 하려면 결과를 지운다.
  - 예: 번들 조정부터 다시 → `bba/`, `images/`, `stereo/`, `final/`을 지운다.
  - 예: 섞는 σ만 바꾸기 → `final/`을 지우고 `STAGES="blend ortho products diag" SIGMA=50 ./run_block.sh ...`
- `undist`와 `stereo`는 파일 단위로 이어서 진행한다. 카메라가 바뀌었다면 `images/`와 `stereo/`를 반드시 지운다.

---

## 5. 결과물

`<제품 폴더>/` (블록 `B`):

| 파일 | 내용 |
|---|---|
| `B_kjr_dsm.tif` | DSM, COG Float32, nodata -9999, 정사영상과 같은 격자 |
| `B_kjr_ortho.tif` | 정사영상, COG RGBA |
| `B_kjr_dsm_envi.tif`, `B_kjr_ortho_envi.tif` | ENVI용 비압축 Float32 (nodata NaN). ENVI는 압축 TIFF를 읽지 못한다. |
| `B_kjr_ortho.jpg`, `B_kjr_ortho_85pct.jpg` | 미리보기 (+ `.wld`) |
| `B_kjr_dsm_hillshade.jpg` | 전체 해상도 힐쉐이드 |
| `B_kjr_dsm_color_hillshade_small.jpg` | 높이 5배 과장 컬러 힐쉐이드 (0.76 m) |

`<작업 폴더>/final/dsm_fused.tif`는 DSM 원본 격자(정사영상 픽셀의 2배, KJ-1은 13 cm)이다.

---

## 6. 결과 확인

`diag` 단계가 로그에 남기는 것:

1. **GCP 높이** (`tools/gcp_check.py`): `ref`, 번들 조정만 한 DSM, 최종 DSM의 GCP 높이 잔차.
   - KJ-1 값: 3.8 / 19.8 / **5.2 cm RMS**.
   - 최종 DSM이 `ref`보다 크게 나쁘면 `SIGMA`를 줄인다.
2. **사진 좌표 무늬** (`tools/diag_imgfixed.py`): `image-locked mean pattern std`와 `even/odd corr`.
   - KJ-1 값: ODM 카메라 1.14 cm / 0.84 → 최종 **0.32 cm / 0.19**.
   - 상관이 0.5를 넘으면 카메라 모델 무늬가 남아 있다는 뜻이다.
3. **번들 조정 로그** (`bba/bba.log`):
   - tie 잔차 중앙값: KJ-1은 0.087 → 0.045 px.
   - GCP 재투영 RMS, 카메라 이동량, 격자 보정 최대값도 기록된다.

수동으로 쓰는 도구:

```bash
source env.sh
$PY tools/diag_consec.py <work>/stereo 53 64     # 한 줄의 연속 사진 쌍 사이 단차 (cm)
$PY tools/color_hillshade.py <dsm.tif> out.jpg --zmin 4.85 --zmax 5.89
$PY tools/gcp_check.py <gcp_list.txt> a.tif b.tif
```

---

## 7. 주요 상수 (코드 안)

| 위치 | 상수 | 값 | 비고 |
|---|---|---|---|
| `kjr/bba.py` | `RSTEP` / `GSTEP` | 30 / 64 px | 반경 / 격자 보정 간격. 격자 30 px는 Schur 행렬이 16 GB를 넘는다 (RAM 30 GB). |
| `kjr/bba.py` | `SIG_TIE` / `SIG_GCP` | 0.1 / 0.3 px | 관측 가중 |
| `kjr/bba.py` | `SIG_ORTH` | 1e-4 | 보정을 저차 모드(반경 1, r, r³…, 격자 3차 이하)와 직교시켜 블록 휨과 섞이지 않게 한다 |
| `kjr/fuse.py` | `ERODE` / `RAMP` | 30 / 60 셀 | 사진 쌍 래스터 가장자리 제거 / 가중 경사 |
| `kjr/fuse.py` | `--sigma` | 1.5 셀 | 최종 Gaussian 평활 |
| `kjr/lowblend.py` | `--sigma` | 30 m | |
| `kjr/ortho.py` | `POW` | 4 | 사진 중심 가중 (1-r)^4 |

---

## 8. 한계와 주의

- 실제로 검증한 블록은 **KJ-1 하나**다 (DJI FC6310, 68장, GCP 5점, 평탄한 사막). 다른 블록에서는 `diag` 결과를 반드시 확인한다.
- KJ-1 번들 조정의 격자 보정 최대값이 15 px까지 커졌다. tie point가 적은 사진 가장자리로 추정한다. DSM에서는 이상을 찾지 못했다.
- KJ-1 최종 DSM 동쪽에 희미한 세로 띠가 하나 남아 있다. 사진 겹침 경계로 추정한다.
- `tools/diag_consec.py`는 DJI 기본 파일명(`DJI_%04d.JPG`)을 가정한다.
- 카메라는 블록당 하나(ODM 재구성의 첫 카메라)라고 가정한다.
