# ODM-KJR 매뉴얼 (기술 문서)

> **사용법만 필요하면 [`manual/README.md`](manual/README.md)** (사용자 매뉴얼, 표준 컨트롤 파일·실행 템플릿은 `manual/templates/`).

드론 사진과 GCP 목록으로 DSM과 정사영상을 만드는 도구다. **ODM은 필요 없다.** 실행은 컨트롤 파일(`*.cfg`)과 `kjr_run`으로 한다 (§4).
사진 위치·자세와 카메라 값(SfM)은 원본 OpenSfM(BSD)으로 구하고(§4.1), 다음 세 단계를 자체 코드로 한다.
끝난 ODM 처리 폴더에서 시작하는 예전 방식도 `SFM_SOURCE`로 쓸 수 있다. 라이선스는 `THIRD_PARTY_LICENSES.md`.

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

`kjr_run`이 아래 단계를 순서대로 실행한다 (`kjr_run --stages`). 번호 1–4는 OpenSfM SfM (§4.1), 5–15는 아래 표, 16 `EXPORT`는 전달 폴더(§5.1).
컨트롤 파일의 `START`, `END`로 범위를 정하고, `EXISTING`으로 이미 있는 결과를 어떻게 할지 정한다 (§4.3).

| 번호 | 단계 | 하는 일 |
|---|---|---|
| 1 | `SFM_FEATURE` | OpenSfM: EXIF, 카메라 모델(EXIF 초점거리에서 시작하는 brown), SIFT 특징점 |
| 2 | `SFM_MATCH` | OpenSfM: 매칭, 트랙 |
| 3 | `SFM_RECON` | OpenSfM: 증분 복원, GCP 포함 번들 조정 |
| 4 | `GEOREF` | ENU → UTM, 정사영상 격자, SfM 점으로 참조 DSM (`osfm_import.py`) |
| 5–15 | `PAIRS` … `DIAG` | 아래 표 |
| 16 | `EXPORT` | 전달 폴더 (`EXPORT_DIR` 필요, 기본 범위 밖) |

SfM을 세 단계로 나눈 이유: 특징점 추출은 설정과 무관하게 재사용할 수 있고, 매칭·트랙은 매칭 설정에만, 복원은 GCP·카메라 설정에만 따른다.
GCP나 카메라 설정만 바꾸면 `START : SFM_RECON`으로 복원만 다시 한다 (JJ-3: 특징점 1분, 매칭 2분, 복원 8분).

| 단계 | 하는 일 | 결과 (`<work>/`) | KJ-3 소요 (RTX 4060 8 GB, 16코어), 최적화 전 → 후 |
|---|---|---|---|
| 5 `PAIRS` | 사진 쌍 목록. `REF_DSM : stereo`이면 SfM 카메라로 stereo + 융합해 긴 파장 참조 DSM도 만든다 | `ref/stereo_meta.json`, `ref/dsm_fused.tif` | `REF_DSM : stereo`일 때 약 45분 → 약 15분, 기본은 수 초 |
| 6 `TIEPOINT` | 원본 JPEG에서 tie point 추출: 사진별 SIFT 1회 → GPU 유도 정합 → MAGSAC++ → GPU ALSC | `ties/<a>__<b>.npz` | 39분 → **8분** (215쌍) |
| 8 `BBA` | 번들 조정 (LM, 점 소거 Schur는 C++/OpenMP 커널, Cauchy 가중) | `bba/bba.json`, `bba/bba.log`, `bba/odm/` | 22분 → **5분** |
| 9 `UNDISTORT` | 새 카메라로 Lanczos 왜곡 보정 | `images/*.JPG.tif` | 1분 |
| 10 `DSM_EXTRACTION` | 사진 쌍별 census-SGM (fsgbm, MGM, CUDA). 탐색 중심은 ODM DSM. CPU 작업은 스레드로 GPU 정합과 겹친다. 서브픽셀 커널과 SGM 실행 루프는 자체 구현(`kjr/fastmatch.py`) | `stereo/pairs/*.npz` | 21분 → 12분 → **6.7분** |
| 11 `FUSE` | 강건 융합 (중앙값 + NMAD, 가장자리·외곽 가중). 셀별 통계는 GPU. 사진 쌍이 1–2개뿐인 셀은 스파이크 기준 12 cm + 평활 쪽으로 혼합 (`lown_smooth`, 끄기 `--no-lown`) | `stereo/dsm_fused.tif` | 17분 → **1분** |
| 12 `BLEND` | 긴 파장(σ 30 m)은 참조 DSM, 세부는 stereo. GCP 대비 평균 높이 편차가 유의하면(> 표준오차 2배) 상수 하나만큼 옮김 | `final/dsm_fused.tif` | 수 초 |
| 13 `ORTHO` | 최종 DSM으로 정사영상 (가장 수직에 가까운 사진에 가중), 타일 병렬 | `final/ortho.tif` | 4분 → **2분** |
| 14 `PRODUCTS` | COG, ENVI, JPG 제품 파일 | `<products>/` | 약 2.5분 (디스크 쓰기) |
| 15 `DIAG` | 컬러 힐쉐이드, GCP 높이, 사진 좌표 무늬 검사 | `diag/`, 로그 | 수십 초 |

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

### 2.2 가장자리 노이즈 (2026-10-27)

사진 쌍이 1–2개뿐인 셀(블록 가장자리, JJ-1 동쪽, JJ-2 북동쪽)은 이상값 검사를 할 수 없어서 미세 노이즈가 다른 곳의 약 2배였다 (고주파 p95 14–17 cm 대 5–6 cm).
`fuse.py`의 `lown_smooth`는 이런 셀을 NaN을 고려한 가우시안(σ 2셀)과 섞는다. 원래 값의 가중치는 쌍 1개면 0.35, 2개면 0.675, 3개 이상이면 1이고, 이 가중치 자체를 다시 평활(σ 8셀)해서 쌍 수가 바뀌는 곳에 경계가 생기지 않게 한다. 처음 판(쌍 1개 = 완전 평활, σ 3)은 JJ-3 서쪽 가장자리에서 흐릿한 사각형 조각이 보여서 바꿨다. 이런 셀은 스파이크 기준도 25 → 12 cm로 낮춘다.
JJ-1 결과(첫 판): 쌍 1개 셀의 고주파 노이즈 p95가 5.4 → 1.5 cm, 쌍 2개 셀은 6.9 → 3.9 cm로 줄었다. 쌍 3개 이상인 셀과 GCP 높이(RMS 8.4 cm)는 그대로다.

### 2.3 다중 사진 tie 트랙 (선택, 2026-10-27)

7 `TRACKS` 단계(`kjr/tracks.py`, `USE_TRACKS : yes`일 때만 실행)는 사진 쌍마다 저장한 SIFT 특징점 번호로 같은 점을 여러 사진에 걸쳐 묶는다. 트랙마다 사진 중심에 가장 가까운 관측을 기준으로 정하고, 나머지 관측은 모두 그 기준에서 ALSC로 다시 잰다.
JJ-3에서 3장 이상 보이는 점이 8 → 42 %로 늘었지만 결과 DSM은 거의 같았다 (GCP 9.9 → 9.7 cm, 사진 좌표 무늬 0.42 → 0.43 cm). 서쪽 가장자리의 가로 물결도 그대로라서, 원인은 번들 조정의 블록 강성이 아니다. KJ-3도 같았다 (stereo GCP 10.5 → 10.2 cm). 시간은 약 10–20분 더 들기 때문에 기본으로는 쓰지 않는다.

### 2.4 GCP 높이 기준 (2026-10-27)

긴 파장은 ODM DSM에서 가져오므로 ODM의 수직 편차(몇 cm)도 함께 들어온다. `lowblend.py --gcp`가 GCP 위치에서 DSM − GCP의 평균을 구하고, 유의할 때만 그만큼 전체를 옮긴다.
6개 블록 leave-one-out 검증 결과:

| 블록 | 보정 없음 | 상수 이동 (LOO) | 평면 (LOO) |
|---|---|---|---|
| KJ-1 | 2.5 | 3.0 | 2.9 |
| KJ-2 | 5.8 | 3.3 | 5.3 |
| KJ-3 | 7.1 | 4.2 | 7.6 |
| JJ-1 | 8.4 | 7.4 | 7.1 |
| JJ-2 | 7.4 | 5.5 | 6.5 |
| JJ-3 | 9.9 | 9.3 | 17.7 |

GCP 5개로 평면을 맞추면 과적합되므로(JJ-3 17.7 cm) 상수 이동만 쓴다. KJ-1은 편차(+0.8 cm)가 유의하지 않아 옮기지 않는다.

---

## 3. 설치

필요한 것:
- Linux, NVIDIA GPU (8 GB면 충분), CUDA 드라이버
- Python 환경 (conda prefix): numpy, scipy, opencv, gdal, matplotlib, pillow, torch, cuda-python
- UAV 정합 라이브러리 `UAV_LIB`: `tiepoints`(SIFT + MAGSAC++ + ALSC)와 `matcher`(census SGM). 현재 위치는 `/work10/planet_stereo/UAV`이다.
- GDAL 명령행 도구: `gdal_translate`, `gdalwarp`, `gdaldem`
- `g++` (OpenMP): 번들 조정의 C++ 커널 `kjr/csrc/schur.cpp`를 처음 실행할 때 `kjr/_schur.so`로 자동 컴파일한다. 직접 하려면 `cd kjr/csrc && g++ -O3 -march=native -fopenmp -shared -fPIC schur.cpp -o ../_schur.so`. 컴파일이 안 되면 scipy 방식으로 돌아간다 (3–4배 느림).
- 원본 OpenSfM (기본 처리 §4.1): `tools/build_opensfm.sh <폴더>`로 별도 conda 환경에 빌드한다. **ODM은 필요 없다.**
- (ODM에서 시작하는 예전 방식만) 끝난 ODM GCP 처리 폴더. 다음 파일이 있어야 한다:
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

블록마다 컨트롤 파일(`*.cfg`, 텍스트) 하나를 쓰고 `kjr_run`으로 실행한다.

```bash
cd /work/eridu_drone/kjr/odm-kjr
./kjr_run cfg/eridu/JJ3.cfg                       # 블록 하나
./kjr_run cfg/eridu/*.cfg                         # 여러 블록을 차례로 (하나가 실패해도 다음 블록 진행)
./kjr_run cfg/eridu/JJ3.cfg --start BBA --end DIAG --existing archive   # 범위·처리 방식을 명령에서 덮어쓰기
./kjr_run cfg/eridu/JJ3.cfg --dry-run             # 검사만: 실행할 단계와 명령을 보여 주고 아무것도 안 함
./kjr_run --stages                                # 단계 목록
./kjr_run --keys                                  # 컨트롤 파일 항목과 기본값
```

오래 걸리므로 `nohup ./kjr_run ... > run.out 2>&1 &`로 띄우고 `<WORK_DIR>/kjr_run.log`를 본다.
**원본 사진과 `SFM_SOURCE`의 ODM 폴더는 읽기만 한다.** 결과는 `WORK_DIR`, `OUTPUT_DIR`, `EXPORT_DIR`에만 쓴다.

### 4.1 컨트롤 파일

```
# JJ3.cfg   '#' 뒤는 주석. 형식 "항목 : 값". 항목 이름은 대소문자 구분 없음
TITLE        : Eridu JJ-3
BLOCK        : JJ3
INCLUDE      : ../common.cfg                     # 공통 설정 (환경 경로 등)
INPUT_DIR    : /work/eridu_drone/JJ-3           # 원본 사진 *.jpg
GCP_FILE     : /work/eridu_drone/final_output_jj3/gcp_list_v2.txt
WORK_DIR     : /work/eridu_drone/kjr/runs_osfm/${BLOCK}
OUTPUT_DIR   : /work/eridu_drone/kjr/products_osfm/${BLOCK}
START        : SFM_FEATURE                       # 단계 이름 또는 번호
END          : DIAG
EXISTING     : skip                              # skip | archive | delete
```

- 상대경로는 그 컨트롤 파일이 있는 폴더 기준이다. `${BLOCK}`처럼 다른 항목 값을 쓸 수 있다.
- `INCLUDE`한 파일의 값은 기본값 역할이다. 같은 항목이 있으면 **포함하는 파일이 이긴다** (줄 위치와 무관).
- 모르는 항목(오타)은 오류로 멈춘다 (`unknown key STRAT (did you mean START?)`). 필수 항목이 없거나, 파일·폴더가 없거나, `END`가 `START`보다 앞이어도 실행 전에 멈춘다.
- 예: `cfg/example.cfg`, `cfg/common.cfg`, 이 현장 6개 블록 `cfg/eridu/*.cfg`.

| 항목 | 기본값 | 의미 |
|---|---|---|
| `BLOCK`, `INPUT_DIR`, `GCP_FILE`, `WORK_DIR`, `OUTPUT_DIR` | (필수) | 블록 이름, 사진 폴더, GCP 목록(OpenSfM/ODM 형식), 작업 폴더, 제품 폴더 |
| `TITLE` | | 로그에 쓰는 설명 |
| `START`, `END` | 첫 단계, `DIAG` | 실행 범위 (이름 또는 번호, §2) |
| `EXISTING` | `skip` | §4.3 |
| `SFM_SOURCE` | `opensfm` | `opensfm` = 1–4단계 실행. 끝난 ODM 폴더를 주면 그 결과에서 시작 (1–4단계 없음, 예전 방식) |
| `OSFM_CONFIG` | `tools/opensfm_config.yaml` | OpenSfM 설정 |
| `OSFM_PROCESSES` | `auto` (CPU − 2) | OpenSfM 프로세스 수 |
| `CAMERA_OVERRIDES` | `auto` | `auto` = 카메라마다 EXIF 초점거리에서 시작하는 brown. 아니면 OpenSfM `camera_models_overrides.json` 파일 |
| `GSD` | `auto` | 정사영상 해상도 (m). `auto` = 비행고도 / 초점거리 |
| `REF_CELL` | 2 | SfM 점으로 만드는 참조 DSM 격자 (m) |
| `REF_DSM` | `auto` | 긴 파장 참조 DSM이자 stereo 탐색 중심. `auto` = SfM DSM (`opensfm`: OpenSfM 점, ODM: ODM DSM). `stereo` = `PAIRS`에서 SfM 카메라 stereo로 만듦 (느림, JJ-1 GCP RMS 13.6 cm 대 `auto` 8.3 cm). 또는 tif 파일 |
| `SHAPE_DSM` | `ref` | `BLEND`의 긴 파장 기준. `ref` = `REF_DSM`, `self` = 긴 파장 보정 안 함 (GCP 오프셋만), 또는 tif |
| `BLEND_SIGMA` | 30 | `BLEND`의 Gaussian σ (m) |
| `BBA_ITERS` | 30 | 번들 조정 반복 횟수 |
| `BBA_OPTIONS` | | `bba.py` 추가 옵션. 예: `--sig-gcp 0.03` |
| `TIE_WORKERS` | 2 | tie point 병렬 작업 수 (GPU 공유, 메모리가 모자라면 자동 재시도) |
| `USE_TRACKS` | `no` | `yes` = 7 `TRACKS` 실행, BBA에 사용 (§2.3) |
| `PRODUCT_NAME` | `${BLOCK}_kjr` | 제품 파일 이름 앞부분: `<이름>_dsm.tif`, `<이름>_ortho.tif`, … |
| `EXPORT_DIR`, `REGIONAL_DEM` | | `EXPORT` 단계 (§5.1) |
| `KJR_ENV`, `UAV_LIB` | `env.sh` 기본 | KJR 파이썬 환경, UAV 정합 라이브러리 |
| `OSFM_ENV`, `OSFM_SRC` | `/work/eridu_drone/kjr/opensfm/{env,src}` | OpenSfM 환경과 소스 (`tools/build_opensfm.sh`) |

### 4.2 SfM: 원본 OpenSfM (1–4단계)

```bash
tools/build_opensfm.sh /work/eridu_drone/kjr/opensfm      # 한 번만: 원본 OpenSfM을 수정 없이 빌드 (65571c8a)
```

| 단계 | 하는 일 | 코드 |
|---|---|---|
| 1–3 | OpenSfM: SIFT, 매칭, 트랙, 증분 복원, GCP 포함 번들 조정. 설정은 ODM 3.5.6이 쓰던 값 (`tools/opensfm_config.yaml`, ODM 전용 dspsift 대신 SIFT) | `kjr/jobctl.py` → `opensfm_main.py` |
| 4 `GEOREF` | OpenSfM 국지 좌표(ENU) → UTM. 점과 카메라 중심은 정확히 변환, 회전은 카메라마다 국지 상사변환. 사진 촬영 범위로 정사영상 격자, OpenSfM 점으로 참조 DSM | `kjr/osfm_import.py` |

- 결과: `<WORK_DIR>/sfm/osfm` (OpenSfM 데이터셋), `<WORK_DIR>/sfm/odm` (5단계 이후의 입력).
- JJ-3 검증 (2026-10-27):
  - OpenSfM 11분, 117장 전부 복원. 초점거리 0.678 (ODM 0.666).
  - 번들 조정 후 GCP 높이 오차가 ODM에서 시작했을 때와 4 cm 이내로 같다.
  - GCP를 하나씩 빼는 검증 (블렌드 전): 102 cm, ODM 시작은 99 cm. 2.2 km 길이에 GCP 5점이라 둘 다 크다.
  - 최종 DSM GCP RMS 5.3 cm (ODM 시작 전달본 9.9 cm). 전달본과의 차이 NMAD 14 cm, 수백 m 파장 (참조 DSM이 ODM 조밀 DSM 대신 OpenSfM 희소점).

### 4.3 범위와 기존 결과 (`START`, `END`, `EXISTING`)

- `START` 앞의 단계는 결과가 있어야 한다. 없으면 "`BBA (before START) has no results … Set START to BBA or earlier`"로 멈춘다.
- `START`–`END` 안의 단계는 `EXISTING`에 따라:

| `EXISTING` | 결과가 이미 있는 단계 |
|---|---|
| `skip` (기본) | 끝난 단계는 건너뛴다. 덜 끝난 단계(`TIEPOINT`, `UNDISTORT`, `DSM_EXTRACTION`)는 남은 것만 이어서 한다 |
| `archive` | 다시 실행한다. 옛 결과는 먼저 `<WORK_DIR>/old/<시각>/`으로 옮긴다 (같은 디스크 안 이동이라 빠르지만 용량은 그대로 차지) |
| `delete` | 다시 실행한다. 옛 결과는 먼저 지운다 |

- **낡은 결과 (stale)**: 어떤 단계를 실행하면 그 뒤 단계 결과는 모두 낡은 것으로 표시된다 (`<WORK_DIR>/job/stale`). 낡은 단계는
  - `skip`이어도 건너뛰지 않는다 (지우고 다시 실행: 이어서 하는 단계가 옛 사진·옛 쌍을 그대로 쓰는 것을 막는다).
  - `START` 앞에 있으면 멈춘다. 예: `BBA`만 다시 돌린 뒤 `START : DSM_EXTRACTION`이면 "`UNDISTORT … is stale`".
  - `END` 뒤에 남으면 끝에 "results of later stages are now stale"로 알려 준다.
- 예:
  - 섞는 σ만 바꾸기: `BLEND_SIGMA : 50`, `START : BLEND`, `END : DIAG`, `EXISTING : archive`.
  - GCP 수정 후: `START : SFM_RECON` (특징점·매칭 재사용).
  - 실패 후 이어서: 같은 파일을 `EXISTING : skip`으로 다시 실행.

### 4.4 기록

실행할 때마다 `<WORK_DIR>/job/<시각>/`에 남긴다:
- `control.cfg` (원본 복사), `resolved.cfg` (기본값까지 채운 모든 값), `git.txt` (커밋, 수정 사항)
- `code/`: 실행한 코드 사본. 모든 단계는 이 사본으로 실행되므로 **실행 중에 저장소를 고쳐도 돌고 있는 작업은 바뀌지 않는다.**

그 밖에 `<WORK_DIR>/kjr_run.log` (단계별 시작·끝·소요 시간·요약), `<WORK_DIR>/job/logs/<번호>_<단계>.log` (단계 전체 출력), `<WORK_DIR>/job/done/` (끝난 단계 표시).

### 4.5 예전 명령 (호환용)

`run.sh <블록> <사진> <GCP> <작업> <제품>`과 `run_block.sh <블록> <ODM 폴더> <사진> <GCP> <작업> <제품>`은 컨트롤 파일을 `<작업>/job/`에 만들어 `kjr_run`을 부른다.
환경변수 `STAGES`(처음–끝 단계를 범위로), `REF_DSM`(`run_block.sh`에서 없으면 `stereo`, 예전 기본과 같음), `SHAPE_DSM`, `SIGMA`, `BBA_ITERS`, `TIE_WORKERS`, `BBA_ARGS`, `GSD`, `EXISTING`, `DRY_RUN=1`을 받는다.
전달본(2026-10-26)은 `backup/delivered_20261026/rerun_delivered.sh`로 다시 만든다 (그때 코드 그대로).

### 4.6 자체 SfM으로 처리 (`run_noodm.sh`) — 실패, 쓰지 말 것

JJ-3에서 초점거리를 10 % 틀리게 추정해 블록 끝에서 수 m 오차가 났다. 기록으로만 남긴다. 대신 `kjr_run`(§4)을 쓴다.

```bash
./run_noodm.sh <블록> <원본 사진 폴더> <GCP 목록> <작업 폴더> <제품 폴더>
```

ODM이 하던 일을 다음 단계로 대신한다.

| 단계 | 하는 일 | 코드 |
|---|---|---|
| `init` | EXIF/XMP에서 초기 카메라: GPS → UTM(경도로 영역 계산), 짐벌 요·피치·롤(없으면 비행 요, 그것도 없으면 GPS 진행 방향), 초점거리(DJI 보정값 → 35 mm 환산 → 초점거리 mm × 센서 해상도). 격자용 빈 정사영상 틀, coords.txt | `kjr/exif_init.py` |
| `pairs` | EXIF 자세로 지상 범위 → 겹치는 사진 쌍 (겹침 5 % 이상) | `stereo.py --meta-only` |
| `ties`, `tracks` | 기존과 같다 (자세와 무관한 정합) | `ties.py`, `tracks.py` |
| `rel` | tie만으로 번들 조정, GPS를 약한 사전정보로 (5 m, 0.1 rad), 왜곡 0 고정 | `bba.py --no-gcp` |
| `abs` | GCP를 삼각측량해 측량 좌표로 유사변환 (높이 기준면 포함), 그다음 GCP와 함께 번들 조정. `CAMERA`가 있으면 내부 파라미터 고정, 없으면 자체 보정 | `bba.py --align-gcp [--free-intr] --sig-gcp 0.03` |
| (선택) `CHECK=1` | GCP를 하나씩 빼고 `abs`를 다시 풀어 검사점 오차 보고 | `bba.py --holdout` |
| `refdsm` | 조정된 tie 점(약 150만 개)으로 2 m 참조 DSM → stereo 탐색 중심. 긴 파장 형태로는 쓰지 않는다 (`SHAPE_DSM=self`: KJ-3 GCP RMS 6.1 → 13.2 cm로 나빠졌음) | `kjr/refdsm.py` |
| 이후 | `run_block.sh`를 `abs/odm`에 대해 실행 (보정 격자 BBA, stereo, 융합, 정사영상, 제품) | |

필요 조건: 카메라 한 종류, 거의 연직인 사진, GPS, 2장 이상에서 보이는 GCP 3개 이상. `RelativeAltitude` 태그가 없으면 `EXIF_ARGS="--agl <비행고도 m>"`.

**정확도 한계 (중요).** 모든 사진이 연직이고 비행선이 평행한 블록에서는 렌즈 왜곡(k1–k3, 초점거리)과 블록 곡률이 서로 맞바뀐다. GCP 5개로는 이 둘을 구분하지 못한다. KJ-3에서 GCP를 하나씩 빼고 검사점으로 쓴 결과:

| 방식 | 검사점 높이 오차 (5점) | RMS |
|---|---|---|
| ODM 없이, 내부 파라미터 자유 자체 보정 (`--free-intr`, GCP 0.03 px) | +0.79, −0.64, +1.17, −0.55, +0.74 m | 약 80 cm |
| 위 + OpenSfM 방식 내부 파라미터 사전정보 (`--intr-prior 0.01,...`) | 같음 (tie가 0.1 px로 매우 정확해서 사전정보가 묻힘) | 약 80 cm |
| **ODM 없이, 같은 기종의 다른 블록(KJ-1) 보정값 고정 (`CAMERA=`)** | +0.43, −0.45, −0.01, −0.36, +0.44 m | **약 37 cm** |
| 참고: ODM 해에서 시작 (같은 블록의 ODM 보정값) | −0.30, −0.10, +0.25, −0.09, −0.23 m | 약 21 cm |

ODM 자신도 블록마다 다른 왜곡 계수를 냈다 (FC6310 다섯 블록: 초점거리 0.6626–0.6685, k2 −0.0045 ~ −0.0086). 같은 모호성을 블록마다 다르게 정한 것이다.
GCP에 큰 가중치를 주면 GCP 자체는 잘 맞지만 (KJ-3 GCP RMS 6.1 cm) GCP 사이 형태는 검사점에서 수십 cm 틀린다. GPS 높이를 카메라 사전정보로 넣는 방법(`--gps-prior`)도 효과가 없었다 (사진별 GPS 높이 잡음 0.6 m).
보정값이 없는 카메라(JJ-3, FC9313: 35 mm 환산값 24 mm로 f 2840 px)는 자체 보정이 초점거리를 10 % 틀리게 잡았다 (0.601 대 ODM 0.666).

그래서:
- 이 블록들처럼 연직만 찍은 자료는 **ODM 해에서 시작하는 `run_block.sh`가 여전히 가장 정확하다.**
- `run_noodm.sh`는 ODM 없이 끝까지 처리할 수 있게 해 준다. 하지만 **카메라 보정값(`CAMERA=`)을 주거나**, 교차·경사 비행선이나 GCP가 많은 블록에서 써야 한다. 결과는 반드시 `CHECK=1`(GCP leave-one-out 검사점)로 확인한다.

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

컨트롤 파일에 `EXPORT_DIR`(와 `REGIONAL_DEM`)을 넣고 `START : EXPORT`, `END : EXPORT`로 실행한다. 블록마다 `EXPORT_DIR/<BLOCK>/`을 만들고, `uav_context_map.png`는 그때까지 `EXPORT_DIR`에 있는 모든 블록으로 다시 그린다.
여러 블록을 한 번에 하는 예전 스크립트:

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
