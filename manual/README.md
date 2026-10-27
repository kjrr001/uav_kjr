# uav_kjr 사용자 매뉴얼

드론 사진과 GCP 목록으로 **DSM과 정사영상**을 만드는 방법을 설명한다. 이 문서만 보고 새 블록을 처리할 수 있도록 썼다.
처리 원리, 개발 기록, 상수는 저장소 맨 위의 [`MANUAL.md`](../MANUAL.md)에 있다.

```
manual/
├── README.md               이 매뉴얼
└── templates/
    ├── block.cfg           표준 컨트롤 파일 (블록 하나, 모든 항목과 설명)
    ├── common.cfg          표준 공통 설정 (환경 경로)
    ├── run_job.sh          표준 실행 스크립트 (배경 실행, 상태, 멈추기)
    └── new_block.sh        새 블록 시작: 위 템플릿으로 프로젝트 폴더와 cfg 를 만든다
```

---

## 1. 한눈에 보기

```bash
R=/work/eridu_drone/kjr/odm-kjr                       # 저장소
$R/manual/templates/new_block.sh ~/proj/site1 BLK1 /data/BLK1/images /data/BLK1/gcp_list.txt "현장1 블록1"
cd ~/proj/site1
vi BLK1.cfg                                          # 필요하면 고친다 (보통 그대로)
./run_job.sh check                                   # 검사만
./run_job.sh                                         # 배경 실행
./run_job.sh status                                  # 진행 상황
ls products/BLK1/                                    # 결과: BLK1_kjr_dsm.tif, BLK1_kjr_ortho.tif ...
```

블록 하나에 약 1–1.5시간 걸린다 (JJ-3, 117장 기준: SfM 12분, 나머지 약 1시간).

---

## 2. 준비

### 2.1 컴퓨터

| 항목 | 조건 |
|---|---|
| OS | Linux |
| GPU | NVIDIA, 메모리 8 GB 이상, CUDA 드라이버 (12.2에서 시험) |
| CPU / RAM | 16코어 / 30 GB에서 시험. 번들 조정에 RAM 수 GB |
| 디스크 | 블록당 작업 폴더 약 10 GB + 제품 약 7 GB. `EXISTING : archive`는 옛 결과를 지우지 않으므로 더 든다 |

### 2.2 소프트웨어 (이 서버에는 설치되어 있다)

| 구성 | 위치 (이 서버) | 새로 설치할 때 |
|---|---|---|
| uav_kjr 저장소 | `/work/eridu_drone/kjr/odm-kjr` | `git clone https://github.com/kjrr001/uav_kjr` |
| KJR 파이썬 환경 | `/work/eridu_drone/kjr/env` | conda: python 3.11, numpy, scipy, opencv, gdal, matplotlib, pillow, torch (cu121), cuda-python. `MANUAL.md` §3 |
| UAV 정합 라이브러리 | `/work10/planet_stereo/UAV` | 복사 |
| OpenSfM (원본, BSD) | `/work/eridu_drone/kjr/opensfm/{env,src}` | `tools/build_opensfm.sh <폴더>` |
| GDAL 명령 | `gdal_translate`, `gdalwarp`, `gdaldem` | 시스템 또는 conda |

경로가 다르면 `common.cfg`의 `KJR_ENV`, `UAV_LIB`, `OSFM_ENV`, `OSFM_SRC`만 고친다.
ODM은 필요 없다 (라이선스: `../THIRD_PARTY_LICENSES.md`).

설치 확인:

```bash
R=/work/eridu_drone/kjr/odm-kjr
$R/kjr_run --stages                                   # 단계 16개가 나오면 됨
/work/eridu_drone/kjr/env/bin/python -c "import torch, cv2; from osgeo import gdal; print(torch.cuda.is_available())"   # True
ls /work/eridu_drone/kjr/opensfm/src/bin/opensfm_main.py
```

### 2.3 입력 자료

**사진**
- 한 폴더에 `*.jpg` (대소문자 무관). 하위 폴더는 보지 않는다.
- 카메라 한 종류. EXIF에 초점거리와 GPS가 있어야 한다 (DJI 기본 사진이면 된다).
- 거의 연직으로 찍은 블록에서 시험했다.

**GCP 목록** (OpenSfM/ODM 형식, 텍스트)

```
EPSG:32638
594954.372 3410119.047 5.651 2407.4 1729.1 DJI_20260713080001_0150_D.JPG p1
594954.372 3410119.047 5.651 2367.3 1178.9 DJI_20260713075946_0149_D.JPG p1
...
```

- 첫 줄: 좌표계 (`EPSG:<번호>` 또는 `WGS84 UTM 38N`).
- 다음 줄부터 관측 하나에 한 줄: `동(E) 북(N) 높이(Z) 사진x 사진y 사진파일이름 [점이름]`.
  - 사진 좌표는 원본 사진의 픽셀 (왼쪽 위가 0, 0).
  - 같은 점을 여러 사진에서 찍었으면 사진마다 한 줄.
- GCP는 3점 이상, 점마다 2장 이상에서 보여야 한다. 블록 양 끝과 가운데에 고루 있을수록 좋다
  (JJ-3: 2.2 km에 5점이면 GCP 사이 오차가 1 m에 이를 수 있다).

---

## 3. 새 블록 처리하기

### 3.1 프로젝트 폴더 만들기

```bash
/work/eridu_drone/kjr/odm-kjr/manual/templates/new_block.sh <프로젝트 폴더> <BLOCK> <사진 폴더> <GCP 파일> ["설명"]
```

- 프로젝트 폴더가 없으면 만든다. `common.cfg`, `run_job.sh`가 없으면 템플릿을 복사한다.
- `<BLOCK>.cfg`를 표준 템플릿으로 만든다. 이미 있으면 덮어쓰지 않는다.
- 사진 수, GCP 줄 수, 좌표계를 보여 주고 cfg를 검사한다.
- 같은 현장의 블록은 같은 프로젝트 폴더에 여러 번 실행해서 넣는다.

결과 폴더 구조:

```
<프로젝트 폴더>/
├── common.cfg        공통 설정
├── BLK1.cfg          블록 컨트롤 파일 (블록마다 하나)
├── BLK2.cfg
├── run_job.sh        실행 스크립트
├── logs/             실행 로그 (run_<시각>.log, latest.log)
├── work/BLK1/        중간 결과
└── products/BLK1/    제품
```

### 3.2 컨트롤 파일 확인

대부분 그대로 쓴다. 보통 고치는 것:

| 상황 | 고칠 항목 |
|---|---|
| 정사영상 해상도를 정하고 싶다 | `GSD : 0.08` (m). 기본은 비행고도/초점거리 |
| 작업·제품 폴더를 다른 디스크에 | `WORK_DIR`, `OUTPUT_DIR` |
| 제품 파일 이름을 바꾸고 싶다 | `PRODUCT_NAME` (기본 `<BLOCK>_kjr`) |
| 전달 폴더까지 만들기 | `END : EXPORT`, `EXPORT_DIR`, (`REGIONAL_DEM`) |

모든 항목은 §5 표와 `templates/block.cfg`의 주석에 있다.

### 3.3 실행

```bash
cd <프로젝트 폴더>
./run_job.sh check      # 검사: 실행할 단계와 명령만 보여 준다
./run_job.sh            # 배경 실행 (터미널을 닫아도 계속)
./run_job.sh status     # 실행 중인지, 마지막 진행 상황, 실패 내용
tail -f logs/latest.log # 진행 계속 보기
./run_job.sh stop       # 멈추기
./run_job.sh fg         # 앞에서 실행 (시험용)
```

- `run_job.sh`는 폴더의 모든 블록 cfg(`common.cfg` 제외)를 차례로 처리한다. 일부만 하려면 맨 위 `CFGS=( BLK1.cfg )`.
- 맨 위 `START`, `END`, `EXISTING`을 채우면 모든 cfg의 값을 이번 실행에서만 덮어쓴다.
- 한 블록이 실패해도 다음 블록으로 넘어가고, 끝에 실패한 블록을 알려 준다.
- **멈춘 뒤나 실패 뒤에는 그냥 다시 `./run_job.sh`** 하면 끝난 단계는 건너뛰고 이어서 한다 (`EXISTING : skip`).
- 저장소 명령을 직접 써도 된다: `kjr_run BLK1.cfg [--start X] [--end Y] [--existing Z] [--dry-run]`.

### 3.4 결과 확인

`./run_job.sh status` 또는 `work/<BLOCK>/kjr_run.log`의 `DIAG` 줄을 본다 (§8).

```
DIAG    .../final/dsm_fused.tif p1:-3.3 p2:-2.8 p3:-5.2 p4:+9.6 p5:-1.3 cm; RMS 5.3 cm, mean -0.6 cm
DIAG    hp DSM std 4.24 cm; image-locked mean pattern std 0.42 cm; even/odd corr 0.26
```

그리고 `products/<BLOCK>/<BLOCK>_kjr_dsm_color_hillshade_small.jpg`를 눈으로 본다 (줄무늬, 동심원, 단차, 가장자리 잡음).

---

## 4. 처리 단계

| 번호 | 단계 | 하는 일 | JJ-3 시간 |
|---|---|---|---|
| 1 | `SFM_FEATURE` | OpenSfM: EXIF, 카메라 모델(EXIF 초점거리에서 시작하는 brown), SIFT 특징점 | 2분 |
| 2 | `SFM_MATCH` | OpenSfM: 사진 간 매칭, 트랙 | 3분 |
| 3 | `SFM_RECON` | OpenSfM: 복원, GCP 포함 번들 조정 (카메라 자체 보정) | 8분 |
| 4 | `GEOREF` | 좌표 변환(ENU → UTM), 정사영상 격자, SfM 점으로 참조 DSM | 10초 |
| 5 | `PAIRS` | 스테레오 사진 쌍 목록 | 10초 |
| 6 | `TIEPOINT` | 원본 사진에서 서브픽셀 tie point (SIFT → MAGSAC++ → ALSC, GPU) | 10–30분 |
| 7 | `TRACKS` | 다중 사진 트랙 (`USE_TRACKS : yes`일 때만) | (20분) |
| 8 | `BBA` | 자체 번들 조정 + 사진 좌표 보정(반경, 2차원 격자) | 5분 |
| 9 | `UNDISTORT` | 새 카메라로 왜곡 보정 사진 (Lanczos) | 1분 |
| 10 | `DSM_EXTRACTION` | 사진 쌍마다 census-SGM DSM (GPU) | 20분 |
| 11 | `FUSE` | 쌍별 DSM 합성 | 2분 |
| 12 | `BLEND` | 긴 파장(> `BLEND_SIGMA`)은 참조 DSM에서, GCP 높이 오프셋 | 10초 |
| 13 | `ORTHO` | 정사영상 | 2분 |
| 14 | `PRODUCTS` | 제품 파일 (COG, ENVI, JPG) | 2분 |
| 15 | `DIAG` | 진단: 컬러 음영기복도, GCP 높이, 사진 좌표 무늬 | 30초 |
| 16 | `EXPORT` | 전달 폴더 (`EXPORT_DIR` 있을 때만, 기본 범위 밖) | 1분 |

- `START`, `END`에는 이름이나 번호를 쓴다 (`START : 8` = `START : BBA`).
- `START` 앞 단계는 결과가 있어야 한다. 없으면 시작 전에 멈추고 어디서 시작할지 알려 준다.
- 단계 목록은 `kjr_run --stages`로도 볼 수 있다.

---

## 5. 컨트롤 파일

### 5.1 쓰는 법

```
# 주석
BLOCK        : BLK1                   # 줄 끝 주석도 된다 ('#' 앞에 공백)
INCLUDE      : common.cfg             # 다른 cfg 를 먼저 읽는다
WORK_DIR     : work/${BLOCK}          # 다른 항목 값 쓰기
```

- 한 줄에 `항목 : 값`. 항목 이름은 영문이고 대소문자는 가리지 않는다.
- 상대경로는 **그 cfg 파일이 있는 폴더** 기준이다.
- `INCLUDE`한 파일의 값은 기본값 역할이다. 같은 항목이 있으면 포함하는 파일이 이긴다 (줄 위치와 무관).
- 값을 비우거나 줄을 지우면 기본값이 된다. 같은 파일에 같은 항목이 두 번 있으면 오류다.
- 모르는 항목은 오류로 멈춘다 (`unknown key STRAT (did you mean START?)`). 오타가 기본값으로 조용히 바뀌는 일이 없다.

### 5.2 항목

| 항목 | 기본값 | 의미 |
|---|---|---|
| `BLOCK` | (필수) | 블록 이름 (영문·숫자) |
| `INPUT_DIR` | (필수) | 사진 폴더 |
| `GCP_FILE` | (필수) | GCP 목록 (§2.3) |
| `WORK_DIR` | (필수) | 작업 폴더 |
| `OUTPUT_DIR` | (필수) | 제품 폴더 |
| `TITLE` | | 로그에 쓰는 설명 |
| `START` / `END` | 첫 단계 / `DIAG` | 실행 범위 |
| `EXISTING` | `skip` | `skip` / `archive` / `delete` (§6) |
| `SFM_SOURCE` | `opensfm` | 끝난 ODM 처리 폴더를 주면 1–4단계 없이 그 결과에서 시작 |
| `CAMERA_OVERRIDES` | `auto` | OpenSfM `camera_models_overrides.json`을 주면 그 카메라 모델로 시작 |
| `OSFM_CONFIG` | `tools/opensfm_config.yaml` | OpenSfM 설정 |
| `OSFM_PROCESSES` | `auto` | OpenSfM 프로세스 수 (CPU − 2) |
| `GSD` | `auto` | 정사영상 해상도 (m) |
| `REF_CELL` | 2 | 참조 DSM 격자 (m) |
| `TIE_WORKERS` | 2 | tie point 병렬 작업 수 |
| `USE_TRACKS` | `no` | 다중 사진 트랙 사용 |
| `BBA_ITERS` | 30 | 번들 조정 반복 횟수 |
| `BBA_OPTIONS` | | `bba.py` 추가 옵션 |
| `REF_DSM` | `auto` | 긴 파장 기준 DSM: `auto` (SfM DSM), `stereo`, 또는 tif |
| `SHAPE_DSM` | `ref` | `BLEND` 긴 파장 형태: `ref`, `self` (보정 안 함), 또는 tif |
| `BLEND_SIGMA` | 30 | 긴 파장 경계 (m) |
| `PRODUCT_NAME` | `${BLOCK}_kjr` | 제품 파일 이름 앞부분 |
| `EXPORT_DIR` | | 전달 폴더 |
| `REGIONAL_DEM` | | 전달 폴더 위치도의 바탕 DEM |
| `KJR_ENV`, `UAV_LIB` | `env.sh` 값 | KJR 환경, 정합 라이브러리 |
| `OSFM_ENV`, `OSFM_SRC` | 이 서버 경로 | OpenSfM 환경, 소스 |

`kjr_run --keys`로 기본값을 확인할 수 있다.

---

## 6. 이미 결과가 있을 때 (`EXISTING`)

| 값 | `START`–`END` 안에서 결과가 이미 있는 단계 |
|---|---|
| `skip` (기본) | 끝난 단계는 건너뛴다. 덜 끝난 단계는 이어서 한다 |
| `archive` | 다시 실행한다. 옛 결과는 `WORK_DIR/old/<시각>/`으로 옮긴다 (나중에 직접 지운다) |
| `delete` | 다시 실행한다. 옛 결과는 지운다 |

**낡은 결과 (stale) 자동 처리.** 어떤 단계를 다시 실행하면 그 뒤 단계 결과는 모두 "낡음"으로 표시된다.
- 낡은 단계는 `skip`이어도 건너뛰지 않고, 지우고 다시 실행한다.
- `START` 앞에 낡은 단계가 있으면 시작 전에 멈춘다. 예: `BBA`만 다시 돌린 뒤 `START : DSM_EXTRACTION` → "`UNDISTORT (before START) is stale … Set START to UNDISTORT or earlier`".
- `END` 뒤에 낡은 단계가 남으면 끝에 "results of later stages are now stale"로 알려 준다.

---

## 7. 자주 하는 작업

| 하려는 일 | cfg 또는 `run_job.sh` 맨 위 설정 |
|---|---|
| 처음부터 끝까지 | `START : SFM_FEATURE`, `END : DIAG` (기본) |
| 실패·중단 뒤 이어서 | 그대로 다시 `./run_job.sh` (`EXISTING : skip`) |
| SfM만 해 보기 (사진·GCP 확인) | `END : GEOREF`. 로그의 `reconstruction: N of M images` 확인 |
| GCP 파일을 고친 뒤 | `GCP_FILE` 고치고 `START : SFM_RECON` (특징점·매칭은 그대로 씀) |
| 긴 파장 경계만 바꾸기 | `BLEND_SIGMA : 50`, `START : BLEND`, `EXISTING : archive` |
| 정사영상 해상도만 바꾸기 | `GSD : 0.05`, `START : GEOREF` (격자가 바뀌므로 뒤 단계 모두 다시) |
| 결과를 남겨 두고 비교 실험 | 새 cfg 를 만들고 `WORK_DIR`, `OUTPUT_DIR`만 다르게. 또는 `EXISTING : archive` |
| 전달 폴더 만들기 | `EXPORT_DIR : /전달/폴더`, `REGIONAL_DEM : dem.tif`, `START : EXPORT`, `END : EXPORT` |
| 끝난 ODM 결과에서 시작 | `SFM_SOURCE : /ODM/처리/폴더`, `START : PAIRS` |
| 디스크 정리 | 끝난 블록의 `work/<BLOCK>/images`, `stereo/pairs`, `ties`, `old/` 를 지워도 제품은 남는다 (다시 돌리려면 그 단계부터) |

여러 블록을 같은 설정으로 다시 할 때는 cfg를 고치지 말고 `run_job.sh` 맨 위 `START=`, `END=`, `EXISTING=`을 쓴다.

---

## 8. 결과물과 판단 기준

### 8.1 제품 (`OUTPUT_DIR`)

| 파일 | 내용 |
|---|---|
| `<이름>_dsm.tif` | DSM. COG Float32, nodata −9999, 정사영상과 같은 격자 |
| `<이름>_ortho.tif` | 정사영상. COG RGBA |
| `<이름>_dsm_envi.tif`, `<이름>_ortho_envi.tif` | ENVI용 비압축 (nodata NaN). ENVI는 압축 TIFF를 못 읽는다 |
| `<이름>_ortho.jpg`, `<이름>_ortho_85pct.jpg` | 미리보기 (+ `.wld`) |
| `<이름>_dsm_hillshade.jpg` | 음영기복도 |
| `<이름>_dsm_color_hillshade_small.jpg` | 높이 5배 과장 컬러 음영기복도 (0.76 m). **먼저 이것을 본다** |

`<이름>` = `PRODUCT_NAME` (기본 `<BLOCK>_kjr`). 좌표계는 사진 위치의 UTM (WGS84). 높이는 GCP의 높이 기준을 따른다.
DSM = DTM으로 본다 (식생이 없는 현장 전제).

전달 폴더 (`EXPORT`): `EXPORT_DIR/<BLOCK>/`에 DSM, 정사영상, 컬러 음영기복도, 같은 크기의 작은 정사영상, 사진 촬영 범위 그림. `EXPORT_DIR/uav_context_map.png`는 그 폴더의 모든 블록 위치도.

### 8.2 판단 기준 (`DIAG`, `kjr_run.log`)

| 항목 | 로그 | 좋음 (JJ-3, KJ-1 값) | 나쁨이면 |
|---|---|---|---|
| SfM 복원 | `GEOREF  reconstruction: 117 of 117 images` | 전부 복원, 두 번째 복원 없음 | 빠진 사진이 많으면 사진 겹침·GCP 확인. `WARNING: second reconstruction` 은 블록이 끊긴 것 |
| 번들 조정 | `BBA  it NN: ... tie \|res\| median 0.08 px ... GCP RMS 1.9 px` | tie 중앙값 0.05–0.1 px | 0.3 px 이상이면 사진 품질·매칭 문제 |
| 최종 GCP 높이 | `final/dsm_fused.tif ... RMS 5.3 cm` | 5–10 cm | 20 cm 이상이면 GCP 파일 확인 (사진 좌표, 점 이름) |
| 사진 좌표 무늬 | `image-locked mean pattern std 0.42 cm; even/odd corr 0.26` | std < 1 cm, corr < 0.5 | corr > 0.5 이면 카메라 모델 무늬가 남음 |

참고: GCP는 처리에 쓰였으므로 GCP 높이 RMS는 독립 검증이 아니다. GCP가 적고 블록이 길면 GCP 사이에서는 오차가 더 클 수 있다.

---

## 9. 로그와 기록

| 위치 | 내용 |
|---|---|
| `<프로젝트>/logs/latest.log` | `run_job.sh` 실행 전체 화면 출력 |
| `WORK_DIR/kjr_run.log` | 단계별 시작·끝·소요 시간·요약 (블록마다, 누적) |
| `WORK_DIR/job/logs/<번호>_<단계>.log` | 단계 전체 출력. 실패 원인은 여기 |
| `WORK_DIR/job/<시각>/` | 실행마다: `control.cfg`(원본), `resolved.cfg`(기본값까지 채운 모든 값), `git.txt`, `code/`(실행한 코드 사본) |
| `WORK_DIR/job/done/`, `job/stale` | 끝난 단계 표시, 낡은 단계 목록 |

실행은 `code/` 사본으로 하므로 실행 중에 저장소를 고쳐도 돌고 있는 작업은 바뀌지 않는다.
나중에 같은 결과를 다시 만들려면 `resolved.cfg`와 `git.txt`의 커밋을 쓴다.

---

## 10. 문제 해결

| 메시지 | 원인과 조치 |
|---|---|
| `unknown key XXX (did you mean ...)` | 항목 이름 오타 |
| `XXX is required` | 필수 항목(`BLOCK`, `INPUT_DIR`, `GCP_FILE`, `WORK_DIR`, `OUTPUT_DIR`) 없음 |
| `INPUT_DIR: no *.jpg` / `GCP_FILE: no such file` | 경로 확인. 상대경로는 cfg 폴더 기준 |
| `END X comes before START Y` | 범위 순서 |
| `X (before START) has no results` | 앞 단계가 안 끝났다. `START`를 X나 그 앞으로 |
| `X (before START) is stale` | 앞 단계를 다시 돌린 뒤라 X가 낡았다. `START`를 X로 |
| `START TRACKS does not apply` | `USE_TRACKS : no`인데 `TRACKS`부터 시작. `EXPORT`는 `EXPORT_DIR` 필요 |
| `OSFM_SRC: no bin/opensfm_main.py` | OpenSfM 설치 경로 (`common.cfg`) |
| `X FAILED (...)` + 마지막 25줄 | `WORK_DIR/job/logs/<번호>_X.log` 전체를 본다 |
| `CUDA out of memory` (TIEPOINT) | 자동 재시도한다. 계속되면 `TIE_WORKERS : 1`. 다른 GPU 작업 확인 (`nvidia-smi`) |
| `SFM_RECON` 뒤 사진이 많이 빠짐 | 겹침이 적은 사진, 다른 카메라 모드 사진이 섞였는지 확인 |
| `이미 실행 중` (`run_job.sh`) | `./run_job.sh status`. 정말 아니면 `logs/run.pid`를 지운다 |

---

## 11. 하지 말 것

- 전달본 백업 `/work/eridu_drone/kjr/backup/delivered_20261026/`은 고치거나 지우지 않는다 (2026-10-26 전달본 재현용, 그 안의 `README_RESTORE.md`).
- 같은 `WORK_DIR`로 두 작업을 동시에 돌리지 않는다.
- GPU 작업 두 개를 동시에 돌리면 둘 다 느려지고 메모리 부족이 날 수 있다.
- `run_noodm.sh`(자체 SfM)는 실패한 방식이다. 쓰지 않는다.
