# Akantu for macOS Apple Silicon

Finite element library for solid mechanics, fracture, and contact simulations on Apple Silicon Macs (M1/M2/M3/M4).

**Original Repository**: https://gitlab.com/akantu/akantu
**License**: GNU Lesser General Public License v3.0 (LGPLv3)
**Copyright**: © 2010-2023 EPFL - LSMS

---

## Installation Options

### Option 1: Pre-built Wheel (Recommended - 10 seconds)

Download and install the pre-compiled wheel:

```bash
# Download
wget https://github.com/allamaprabhuani/akantu-macos-apple-silicon/releases/download/v5.0.7-macos-wheel/akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl

# Install Homebrew dependencies
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis

# Install wheel
pip install --user akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl

# Verify
python3 -c "import akantu; print('Installation successful')"
```

**Requirements**: macOS 11.0+, Apple Silicon, Python 3.10

**Documentation**: See [docs/QUICK_START_WHEEL.md](docs/QUICK_START_WHEEL.md)

### Option 2: Build from Source (1-2 hours)

For custom configurations or other Python versions.

#### Step 1: Install Dependencies

```bash
brew install gcc boost scotch eigen open-mpi openblas scalapack cmake
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis
```

If MUMPS fails, edit the formula to use MacPorts mirror:
```bash
brew edit brewsci-mumps
# Change URL to: https://distfiles.macports.org/mumps/MUMPS_5.6.2.tar.gz
# Change SHA256 to: 13a2c1aff2bd1aa92fe84b7b35d88f43434019963ca09ef7e8c90821a8f1d59a
```

#### Step 2: Clone and Build

```bash
git clone https://github.com/allamaprabhuani/akantu-macos-apple-silicon.git
cd akantu-macos-apple-silicon
mkdir build && cd build

CC=gcc-15 CXX=g++-15 FC=gfortran-15 MUMPS_DIR=/opt/homebrew/opt/brewsci-mumps cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_PYTHON_INTERFACE=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=ON \
  -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"

make -j$(sysctl -n hw.ncpu)
sudo make install
```

> **Important:** `-DCMAKE_BUILD_TYPE=Release` is required. Without
> it CMake leaves the build type empty, which means no `-O3`,
> no `-DNDEBUG`, no `-DAKANTU_NDEBUG`, and every internal assertion
> fires in the hot element-wise loops. Measured effect: roughly a
> 20-30x slowdown on phase-field fracture benchmarks.
>
> Example on the Kalthoff-Winkler benchmark (29k nodes, 1500
> explicit steps, Apple M4 Pro): debug build 3552 ms/step, Release
> build 156 ms/step - a 22.8x speedup for a one-flag change.

#### Step 3: Fix Library Paths

```bash
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib
```

#### Step 4: Setup Python for conda base

The compiled `.so` is installed to `/usr/local/lib/python3.10/site-packages/akantu/`. Copy it directly into conda's site-packages so no `PYTHONPATH` override is needed:

```bash
cp /usr/local/lib/python3.10/site-packages/akantu/py11_akantu.cpython-310-darwin.so \
   ~/miniconda3/lib/python3.10/site-packages/akantu/
```

Then create conda activate/deactivate scripts so the library paths are set only when conda base is active (avoids polluting other environments):

```bash
mkdir -p ~/miniconda3/etc/conda/activate.d
mkdir -p ~/miniconda3/etc/conda/deactivate.d

cat > ~/miniconda3/etc/conda/activate.d/akantu_libs.sh << 'EOF'
#!/bin/bash
# libscotch on macOS uses flat namespace and requires libscotcherr to be preloaded
export DYLD_INSERT_LIBRARIES="/opt/homebrew/opt/scotch/lib/libscotcherr.dylib:/opt/homebrew/opt/scotch/lib/libscotcherrexit.dylib"
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:/usr/local/lib${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"
EOF

cat > ~/miniconda3/etc/conda/deactivate.d/akantu_libs.sh << 'EOF'
#!/bin/bash
unset DYLD_INSERT_LIBRARIES
unset DYLD_LIBRARY_PATH
EOF

chmod +x ~/miniconda3/etc/conda/activate.d/akantu_libs.sh
chmod +x ~/miniconda3/etc/conda/deactivate.d/akantu_libs.sh
```

Reload conda and verify:

```bash
conda activate base
python3 -c "import akantu; print('Version:', akantu.__version__); print('MPI:', akantu.has_mpi())"
```

> **Note**: Do not use `conda run` to invoke Python with akantu — macOS strips `DYLD_*` variables from subprocesses launched by `conda run`. Always activate the environment in your terminal first.

> **Note**: Do not add `DYLD_INSERT_LIBRARIES` or `DYLD_LIBRARY_PATH` to `~/.zshrc` globally — this affects all processes and other conda environments. The `activate.d` script scopes them to conda base only.

---

## Features

- **Parallel Computing**: MPI support via Open-MPI
- **Sparse Solvers**: MUMPS direct solver integration
- **Mesh Partitioning**: Scotch library for domain decomposition
- **Python Bindings**: Full Python interface with NumPy/SciPy
- **Phase Field Fracture**: Enhanced models with tension-compression asymmetry for quasi-brittle materials

## Usage Example

```python
import akantu
import numpy as np

# Create mesh
mesh = akantu.Mesh(spatial_dimension=2)
mesh.read("mesh_file.msh")

# Initialize model
model = akantu.SolidMechanicsModel(mesh)
model.initFull()

# Set material properties
mat = model.getMaterial(0)
mat.setParam("E", 210e9)   # Young's modulus
mat.setParam("nu", 0.3)     # Poisson's ratio

# Run simulation
model.solveStep()
```

## Create Custom Wheel

Build a wheel for your Python version:

```bash
python3 packaging/create_standalone_wheel.py
```

Output in `dist/akantu-5.0.7-cp3XX-none-macosx_26_0_arm64.whl`

See [docs/HOW_TO_RELEASE.md](docs/HOW_TO_RELEASE.md) for details.

## Documentation

- **Quick Wheel Installation**: [docs/QUICK_START_WHEEL.md](docs/QUICK_START_WHEEL.md)
- **Wheel Technical Details**: [docs/README_WHEEL.md](docs/README_WHEEL.md)
- **Official Akantu Docs**: https://akantu.readthedocs.io

## Troubleshooting

### Wheel Installation

**ImportError: symbol not found**
- Ensure Homebrew dependencies are installed
- Check Python version matches wheel (3.10)

**Wrong Python version**
```bash
python3 --version  # Check version
# If not 3.10, build custom wheel with packaging/create_standalone_wheel.py
```

### Source Build

**MUMPS not found**
- Set `MUMPS_DIR=/opt/homebrew/opt/brewsci-mumps` in cmake command
- Verify: `ls /opt/homebrew/opt/brewsci-mumps/lib/`

**Scotch errors at runtime** (`symbol not found in flat namespace '_SCOTCH_errorPrint'`)
- `libscotch` on macOS does not embed its dependency on `libscotcherr` — it must be preloaded
- Ensure the conda `activate.d` script is in place (Step 4) and that you ran `conda activate base` in your terminal
- Do **not** use `conda run` — macOS strips `DYLD_*` variables from subprocesses launched that way
- Verify: `echo $DYLD_INSERT_LIBRARIES` (should show the scotcherr paths after activation)

**Python import fails**
- Verify the `.so` was copied: `ls ~/miniconda3/lib/python3.10/site-packages/akantu/py11_akantu*.so`
- Verify conda activation scripts exist: `ls ~/miniconda3/etc/conda/activate.d/`

## Attribution & Citation

This fork enables Akantu on macOS Apple Silicon with phase field enhancements. Original work by EPFL LSMS.

If using in research, please cite:

```bibtex
@article{akantu2024,
  title={Akantu: an HPC finite-element library for contact and dynamic fracture simulations},
  author={Richart, Nicolas and Molinari, Jean-François and others},
  journal={The Journal of Open Source Software},
  year={2024},
  publisher={The Open Journal}
}
```

## System Requirements

- **OS**: macOS 11.0+ (Big Sur or later)
- **CPU**: Apple Silicon (M1/M2/M3/M4)
- **Python**: 3.10 (wheel) or 3.8+ (source build)
- **Disk**: ~5 GB for build, ~500 MB installed

## License

GNU Lesser General Public License v3.0 (LGPLv3)

See [COPYING](COPYING) and [COPYING.lesser](COPYING.lesser) for details.

## Links

- **This Fork**: https://github.com/allamaprabhuani/akantu-macos-apple-silicon
- **Original**: https://gitlab.com/akantu/akantu
- **Documentation**: https://akantu.readthedocs.io
