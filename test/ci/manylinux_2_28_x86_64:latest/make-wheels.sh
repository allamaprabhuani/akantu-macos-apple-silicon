#!/usr/bin/env bash
set -eo pipefail

set +x

export PLAT=manylinux_2_28_x86_64

source /etc/profile

function repair_wheel {
  wheel="$1"
  if ! auditwheel show "$wheel"; then
    echo "Skipping non-platform wheel $wheel"
  else
    auditwheel repair "$wheel" --plat "$PLAT"
  fi
}

# Compile wheels
for PYVER in cp310-cp310 cp311-cp311 cp312-cp312 cp313-cp313 cp39-cp39; do
  PYBIN=/opt/python/${PYVER}/bin
  ccache --zero-stats
  echo "${PYBIN}/pip" wheel . --no-deps -w dist/
  "${PYBIN}/pip" wheel . --no-deps -w dist/
  ccache --show-stats
done

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CI_AKANTU_INSTALL_PREFIX}/lib64

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
  echo repair_wheel "$whl"
  repair_wheel "$whl"
done
