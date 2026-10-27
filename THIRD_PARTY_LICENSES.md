# 외부 소프트웨어와 라이선스 (2026-10-27 확인)

기본 처리(`kjr_run`)는 **ODM을 쓰지 않는다.** ODM(AGPL-3.0) 코드를 복사하거나 import한 곳이 없고,
ODM 프로그램을 실행하지도 않는다. `odm_georeferencing/`, `odm_orthophoto/` 같은 폴더 이름은
파일 배치 형식만 따른 것이며, 내용은 `kjr/osfm_import.py`가 만든다.

| 구성 요소 | 쓰는 곳 | 라이선스 | 배포할 때 할 일 |
|---|---|---|---|
| 이 저장소 (kjr/, tools/, run*.sh) | 전부 | 미정 (LICENSE 파일 없음) | 라이선스를 정해 LICENSE 추가 |
| OpenSfM (mapillary/OpenSfM, 65571c8a, 수정 없음) | SfM (`kjr_run` 1–3단계, `kjr/jobctl.py`) | BSD-2-Clause | 저작권·라이선스 문구 유지. 우리 저장소에는 코드를 넣지 않고 `tools/build_opensfm.sh`가 내려받는다 |
| └ OpenSfM 내장: akaze, vlfeat, pybind11 | OpenSfM 빌드 | BSD 계열 | 위와 같음 |
| └ OpenSfM 내장: openMVS `Interface.h` | `export_openmvs` 명령에만 (우리는 안 씀) | OpenSfM 표기 BSD | – |
| └ OpenSfM 의존: fpdf2 | 보고서 PDF (우리는 안 씀) | LGPL-3.0 | 수정 없이 쓰면 문제 없음 |
| Ceres Solver 2.1 | OpenSfM 번들 조정 | BSD-3 | – |
| SuiteSparse (conda-forge) | Ceres의 희소 행렬 | 모듈별 LGPL/GPL (CHOLMOD 일부 GPL) | 내부 사용은 무관. **빌드된 OpenSfM 바이너리를 남에게 배포**하면 GPL 의무 검토, 또는 SuiteSparse 없이 Ceres 빌드 |
| UAV 정합 라이브러리 (`UAV_LIB`) | tie point, census-SGM | 자체 코드 (LICENSE 없음) | 자체 결정 |
| └ third_party RAFT-Stereo, IGEV | 이 처리에서는 안 씀 | MIT | – |
| numpy, scipy, PyTorch | 전반 | BSD | – |
| OpenCV | 전반 | Apache-2.0 | – |
| GDAL | 래스터 입출력 | MIT | – |
| cuda-python, CUDA NVRTC | GPU 커널 | NVIDIA 라이선스 (무료 사용·재배포 조건 있음) | 재배포 시 NVIDIA 조건 확인 |

## ODM (AGPL-3.0)은 이제 어디에 남아 있나
- **기본 처리에는 필요 없다.** `kjr_run` (컨트롤 파일, 기본 `SFM_SOURCE : opensfm`) = OpenSfM → KJR 단계.
- `SFM_SOURCE : <끝난 ODM 폴더>`로 ODM 결과에서 시작할 수도 있다 (예전 방식). ODM의 **결과 파일**만 읽으므로
  AGPL은 우리 코드에 적용되지 않는다.
- 2026-10-26 전달본은 ODM 결과에서 시작했다 (`backup/delivered_20261026`, 태그 `delivered-20261026`).
  결과물(DSM, 정사영상)에는 라이선스 의무가 없다.

법률 자문이 아니다. 판매·계약 단계에서는 전문가 확인을 받는다.
