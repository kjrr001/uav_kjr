#!/usr/bin/env bash
# Build upstream OpenSfM (mapillary/OpenSfM, BSD-2) into its own conda env, unmodified.
# Usage: build_opensfm.sh <dir>   -> <dir>/src (checkout), <dir>/env (conda prefix)
#   then: export OSFM_ENV=<dir>/env OSFM_SRC=<dir>/src
# Pinned to the commit tested on JJ-3 (2026-10-27). Eigen 3.4.0 (Ceres 2.1 config requires it),
# OpenCV 4 (OpenSfM's CMake looks for calib3d, gone in OpenCV 5).
set -euo pipefail
[ $# -eq 1 ] || { sed -n 2,6p "$0"; exit 1; }
D=$(realpath -m "$1"); COMMIT=65571c8a
mkdir -p "$D"
[ -d "$D/src" ] || git clone --recursive https://github.com/mapillary/OpenSfM "$D/src"
git -C "$D/src" checkout -q $COMMIT && git -C "$D/src" submodule update -q --init --recursive
[ -d "$D/env" ] || mamba create -y -q -p "$D/env" -c conda-forge python=3.10 cmake=3.31 make \
    "libopencv<5=*headless*" "py-opencv<5" ceres-solver=2.1 eigen=3.4.0 llvm-openmp cxx-compiler gflags glog \
    suitesparse pybind11 scikit-build-core numpy scipy pyproj networkx exifread fpdf2 joblib matplotlib \
    pillow python-dateutil pyyaml xmltodict cloudpickle flask
cd "$D/src"
PATH=$D/env/bin:$PATH CMAKE_PREFIX_PATH=$D/env "$D/env/bin/pip" install --no-build-isolation --no-deps .
"$D/env/bin/python" -c "import opensfm.pybundle; print('OpenSfM built:', '$D')"
