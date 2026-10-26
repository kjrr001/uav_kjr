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
| 번들 조정 후 블록 전체가 휨 (±30 cm) | GCP 5점에 비해 보정 자유도가 크다 | 30 m보다 긴 파장은 참조 DSM(권장: ODM 자체 DSM)에서 가져온다 (`lowblend`) |
| 지형이 높거나 낮은 곳이 노이즈로 덮임 (JJ-1 북쪽) | 정합 탐색 범위가 ±8 px뿐이다. 높이로는 ±1 m(기선 130 m), 기선이 긴 쌍은 ±0.5 m다. 처음에는 블록 전체를 평면 하나로 보고 그 주위만 찾았다. | 탐색 중심을 ODM DSM(1 m 평활)에 맞춘다 (`stereo.py --prior-dsm`). JJ-1 북쪽 쌍에서 ODM과 50 cm 넘게 다른 셀이 36% → 2.7% |

시도했다가 버린 방법:
- **내부표정(f, k1–k3)을 풀어서 추정**: 1 m 크기의 돔이 생겼다. 그래서 내부표정은 ODM 값으로 고정한다.
- **사진 쌍별 평면 보정, 쌍 전체 동시 보정**: 오히려 단차를 만들었다. 그래서 `fuse.py`의 기본값을 `--iter 0`으로 둔다.

---

## 2. 처리 단계

`run_block.sh` 하나가 아래 단계를 순서대로 실행한다. 단계마다 결과가 이미 있으면 건너뛴다. 그래서 중간에 실패해도 다시 실행하면 이어서 진행한다.

| 단계 | 하는 일 | 결과 (`<work>/`) | KJ-3 소요 (RTX 4060 8 GB, 16코어), 최적화 전 → 후 |
|---|---|---|---|
| `ref` | ODM 카메라로 stereo + 융합 → 사진 쌍 목록, 긴 파장 참조 DSM. `REF_DSM`을 주면 사진 쌍 목록만 | `ref/stereo_meta.json`, `ref/dsm_fused.tif` | 약 45분 → 약 15분 (`REF_DSM`이면 수 초) |
| `ties` | 원본 JPEG에서 tie point 추출: 사진별 SIFT 1회 → GPU 유도 정합 → MAGSAC++ → GPU ALSC | `ties/<a>__<b>.npz` | 39분 → **8분** (215쌍) |
| `bba` | 번들 조정 (LM, 점 소거 Schur는 C++/OpenMP 커널, Cauchy 가중) | `bba/bba.json`, `bba/bba.log`, `bba/odm/` | 22분 → **5분** |
| `undist` | 새 카메라로 Lanczos 왜곡 보정 | `images/*.JPG.tif` | 1분 |
| `stereo` | 사진 쌍별 census-SGM (fsgbm, MGM, CUDA). 탐색 중심은 ODM DSM. CPU 작업은 스레드로 GPU 정합과 겹친다. 서브픽셀 커널과 SGM 실행 루프는 자체 구현(`kjr/fastmatch.py`) | `stereo/pairs/*.npz` | 21분 → 12분 → **6.7분** |
| `fuse` | 강건 융합 (중앙값 + NMAD, 가장자리·외곽 가중). 셀별 통계는 GPU | `stereo/dsm_fused.tif` | 17분 → **1분** |
| `blend` | 긴 파장(σ 30 m)은 `ref`, 세부는 `stereo` | `final/dsm_fused.tif` | 수 초 |
| `ortho` | 최종 DSM으로 정사영상 (가장 수직에 가까운 사진에 가중), 타일 병렬 | `final/ortho.tif` | 4분 → **2분** |
| `products` | COG, ENVI, JPG 제품 파일 | `<products>/` | 약 2.5분 (디스크 쓰기) |
| `diag` | 컬러 힐쉐이드, GCP 높이, 사진 좌표 무늬 검사 | `diag/`, 로그 | 수십 초 |

블록 하나(`REF_DSM` 사용)에 KJ-3 기준 **약 1시간 47분 → 32분 → 약 27분** 걸린다. 결과 품질은 같다 (KJ-3 GCP 높이 RMS 7.1 → 7.1 cm, 사진 좌표 무늬 0.51 cm / 상관 0.27 → 같음). JJ-3처럼 사진 쌍이 많은 블록(961쌍)은 stereo가 가장 길다 (약 45분 → 약 25분 예상).
디스크는 작업 폴더에 약 10 GB(+ tie 단계 동안 SIFT 캐시 사진당 약 80 MB, 끝나면 지움), 제품 파일에 약 11 GB가 든다.

### 2.1 속도 최적화 기록 (2026-10-26)

| 단계 | 병목 (프로파일) | 바꾼 것 | 결과 동일성 |
|---|---|---|---|
| `fuse` | 셀별 `np.sort` 76 %, 7×7 중앙값 필터·FillNodata를 두 번 계산 | 사진 쌍 전처리 병렬 1회 + 캐시, 중앙값/MAD/가중평균을 GPU(torch), 7×7 중앙값 GPU, 마무리 1회, 형태 연산을 cv2·레이블링으로 | 최대 차이 0.003 mm |
| `bba` | scipy 희소 행렬 곱으로 Schur 보완 (반복당 58초 중 40초), 수치 미분 17회 투영 | `kjr/csrc/schur.cpp`: 점마다 Schur 보완을 직접 누적 (Ceres 방식, OpenMP), 연쇄 법칙 미분 (투영 3회), 고정된 내부 파라미터 열 생략 | 카메라 위치 차이 0.2 mm |
| `ties` | FLANN 근사 최근접 (사진 쌍당 21–29초, 77 %), 사진마다 SIFT 반복 | 사진별 SIFT 캐시, 강한 8000점(안 되면 30000점)으로 호모그래피 → 8×6 타일별 예측 위치 ±120 px 안에서 GPU 정확 2-NN | tie 수 약 1.8배, 최종 DSM 차이 NMAD 0.44 cm |
| `stereo` | CPU 준비/삼각측량이 GPU 정합과 직렬 | 사진 쌍 3개를 스레드로, GPU 정합만 잠금으로 하나씩, 회색조 영상 캐시 | 비트 단위 동일 |
| `ortho` | 타일 직렬 | 타일 스레드 병렬 | 비트 단위 동일 |
| `stereo` (2차) | GPU 정합 쌍당 3.9초 중 서브픽셀 보정 커널이 GPU 시간의 79 % (점마다 쌍선형 보간 7,200회, 시그모이드 13,000회), SGM이 작은 커널 67,890회를 파이썬에서 실행 (실행 오버헤드 > GPU 일) | `kjr/fastmatch.py`: 왼쪽 descriptor는 변위와 무관 → 화소당 1회 계산한 16장 지도; G/H를 안 쓰는 평가는 ±h 샘플 생략; 정수 행이라 쌍선형의 아래 행 생략; 배열 없이 바로 누적 (레지스터 124 → 48개). SGM 단계 루프를 C++(`kjr/csrc/sgmloop.cpp`)로 | 같은 IEEE 연산·순서, 쌍 단독 실행에서 비트 단위 동일 (쌍당 3.9 → 2.1초) |

번들 조정을 일찍 멈추는 것(상대 비용 변화 3e-4)도 시험했으나 카메라 위치가 1.3 cm, 보정 격자가 0.2 px 달라져서 쓰지 않는다 (`--tol` 기본 1e-5 유지).
이전 방식은 옵션으로 남겼다: `bba.py --sparse-schur`, `ties.py --flann`, `fuse.py --cpu`, `stereo.py --uav-kernels`.
빠른 수학 함수(`__expf`, `__fdividef`)로 서브픽셀을 0.3초 더 줄일 수 있지만 결과가 비트 단위로 같지 않아 쓰지 않는다.

주의 (원래부터 있던 비결정성): fsgbm은 실행 순간의 GPU 여유 메모리로 띠(band) 크기를 정한다. 그래서 같은 쌍이라도 다른 쌍과 함께 돌 때와 혼자 돌 때 띠 경계가 달라지고, 몇 개 쌍(KJ-3 215쌍 중 3개)에서 일부 셀이 최대 15 cm 다를 수 있다. 융합의 중앙값이 이를 흡수한다.

---

## 3. 설치

필요한 것:
- Linux, NVIDIA GPU (8 GB면 충분), CUDA 드라이버
- Python 환경 (conda prefix): numpy, scipy, opencv, gdal, matplotlib, pillow, torch, cuda-python
- UAV 정합 라이브러리 `UAV_LIB`: `tiepoints`(SIFT + MAGSAC++ + ALSC)와 `matcher`(census SGM). 현재 위치는 `/work10/planet_stereo/UAV`이다.
- GDAL 명령행 도구: `gdal_translate`, `gdalwarp`, `gdaldem`
- `g++` (OpenMP): 번들 조정의 C++ 커널 `kjr/csrc/schur.cpp`를 처음 실행할 때 `kjr/_schur.so`로 자동 컴파일한다. 직접 하려면 `cd kjr/csrc && g++ -O3 -march=native -fopenmp -shared -fPIC schur.cpp -o ../_schur.so`. 컴파일이 안 되면 scipy 방식으로 돌아간다 (3–4배 느림).
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

권장 실행 (ODM DSM을 참조로):

```bash
REF_DSM=/work/eridu_drone/odm_data/KJ1_gcp/odm_dem/dsm.tif ./run_block.sh KJ1 ...
```

오래 걸리므로 `nohup ... &`로 띄우고 `<작업 폴더>/run_block.log`를 본다.

**ODM 폴더와 원본 사진은 읽기만 한다.** 모든 결과는 작업 폴더와 제품 폴더에만 쓴다.

### 환경변수 옵션

| 변수 | 기본값 | 의미 |
|---|---|---|
| `STAGES` | 전체 | 실행할 단계만 지정. 예: `STAGES="blend ortho products"` |
| `REF_DSM` | (없음) | 긴 파장 참조 DSM을 직접 지정한다. **권장값: ODM 자체 DSM (`<ODM 처리 폴더>/odm_dem/dsm.tif`).** 그러면 `ref` 단계는 사진 쌍 목록만 만들고 약 45분을 아낀다. JJ-1에서 GCP 높이 RMS는 ODM 카메라로 만든 참조를 쓸 때 13.6 cm, ODM DSM을 쓸 때 8.3 cm였고, 무늬 제거 효과는 같았다 (0.46/0.18 → 0.47/0.19). |
| `SIGMA` | 30 | `blend`의 Gaussian σ (m) |
| `BBA_ITERS` | 30 | 번들 조정 반복 횟수 |
| `TIE_WORKERS` | 2 | tie point 사진 쌍 병렬 작업 수(GPU 공유). GPU 메모리가 모자라면 자동으로 재시도한다. |

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

### 5.1 전달용 폴더 만들기

```bash
tools/export_products.sh <products 루트> <작업 루트> <전달 폴더> <지역 DEM.tif> JJ1 JJ2 JJ3 KJ1 KJ2 KJ3
```

블록마다 `<B>_kjr_dsm.tif`, `<B>_kjr_ortho.tif`, `<B>_kjr_dsm_color_hillshade_small.jpg`, `<B>_kjr_ortho_small.jpg`(힐쉐이드와 같은 크기), `<B>_image_footprints.png`(작은 정사영상 위에 블록 경계와 사진별 지상 촬영 범위)를 복사·생성하고, 전달 폴더에 `uav_context_map.png`(지역 DEM 컬러 힐쉐이드 위에 블록 경계와 이름, 경위도 축; 왼쪽 전체, 오른쪽 확대)를 만든다.

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
