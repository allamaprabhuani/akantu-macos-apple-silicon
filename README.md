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
  -DPYTHON_EXECUTABLE=$(which python) \
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

#### Step 4: Setup Python (if enabled)

##### Fix Akantu Python Module Import

If `import akantu` fails due to the module structure, you may need to manually copy the `__init__.py` file so Python can correctly detect the package.

This depends on whether you are using **system Python** or **Miniconda Python**.

---

##### a. Standard Python (System / Homebrew)

If your Python interpreter is located in `/usr/local/bin/python3`.

###### Fix module structure

```bash
sudo cp /usr/local/lib/python3.10/site-packages/akantu/akantu/__init__.py \
        /usr/local/lib/python3.10/site-packages/akantu/__init__.py
```

###### Add to shell profile (`~/.zshrc`)

```bash
export DYLD_INSERT_LIBRARIES="/opt/homebrew/lib/libscotcherr.dylib:/opt/homebrew/lib/libscotcherrexit.dylib"
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:/usr/local/lib:${DYLD_LIBRARY_PATH}"
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:${PYTHONPATH}"
```

---

##### b. Miniconda Python

If your Python interpreter is located in `~/miniconda3/bin/python`.

###### Fix module structure

```bash
cp $HOME/miniconda3/lib/python3.10/site-packages/akantu/akantu/__init__.py \
   $HOME/miniconda3/lib/python3.10/site-packages/akantu/__init__.py
```

###### Add to shell profile (`~/.zshrc`)

```bash
export DYLD_INSERT_LIBRARIES="/opt/homebrew/lib/libscotcherr.dylib:/opt/homebrew/lib/libscotcherrexit.dylib"
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:/usr/local/lib:${DYLD_LIBRARY_PATH}"
export PYTHONPATH="$HOME/miniconda3/lib/python3.10/site-packages:${PYTHONPATH}"
```

---

##### Reload shell

```bash
source ~/.zshrc
```

---

##### Verify installation

```bash
python3 -c "import akantu; print('Success')"
```

---

##### Check which Python you are using

```bash
which python
```

Typical outputs:

```
/usr/local/bin/python3
```

or

```
/Users/username/miniconda3/bin/python
```
```

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

**Scotch errors at runtime**
- Set environment variables in Step 4 above
- Verify: `echo $DYLD_INSERT_LIBRARIES`

**Python import fails**
- Check PYTHONPATH: `echo $PYTHONPATH`
- Verify module location: `ls /usr/local/lib/python3.10/site-packages/akantu/`

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
