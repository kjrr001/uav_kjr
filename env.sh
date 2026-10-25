# source this. Paths can be overridden from the environment before sourcing.
#   KJR_ENV  python env (conda prefix) with numpy scipy opencv gdal torch cuda-python
#   UAV_LIB  UAV matcher library (tiepoints: SIFT+MAGSAC+ALSC, matcher: census SGM)
# Driver CUDA 12.2 -> torch cu121 + cuda-python 12.6, NVRTC 12.1 from the torch
# wheel (a newer NVRTC emits PTX the 535 driver rejects).
export KJR_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KJR_ENV=${KJR_ENV:-/work/eridu_drone/kjr/env}
export UAV_LIB=${UAV_LIB:-/work10/planet_stereo/UAV}
export PY=$KJR_ENV/bin/python
export LD_LIBRARY_PATH=$KJR_ENV/lib/python3.11/site-packages/nvidia/cuda_nvrtc/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=$UAV_LIB:$KJR_REPO
export PROJ_DATA=$KJR_ENV/share/proj   # without it this GDAL opens the dir itself as proj.db
