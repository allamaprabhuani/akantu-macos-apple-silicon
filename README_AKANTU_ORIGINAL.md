# Akantu: macOS Apple Silicon Build

This repository is a fork of the **Akantu** finite element library (Swiss-Made Open-Source Finite-Element Library), cloned from the `features/52-anisotropic-and-at1-phase-field-models` branch at commit `8c5e82d86`.

**Original Repository**: https://gitlab.com/akantu/akantu
**Original Branch**: `features/52-anisotropic-and-at1-phase-field-models`
**Original Commit**: 8c5e82d86

This fork contains platform-specific modifications required for building on macOS with Apple Silicon (arm64), along with enhancements to the phase field fracture framework for modeling tension-compression asymmetry in quasi-brittle materials.

## Building on macOS Apple Silicon

### System Requirements

- **macOS**: 26.0.1 or later (Apple Silicon)
- **Architecture**: arm64
- **Compiler**: GCC 15.2.0 (Homebrew)
- **Python**: 3.10 or later

### Dependencies Installation

Install all required dependencies via Homebrew:

```bash
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis
```

### Build Instructions

1. **Clone this repository**:
```bash
git clone https://github.com/allamaprabhuani/akantu-macos-apple-silicon
cd akantu-macos-apple-silicon
```

2. **Configure CMake**:
```bash
mkdir build && cd build

cmake .. \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=ON \
  -DMUMPS_INCLUDE_DIR=/opt/homebrew/opt/brewsci-mumps/include \
  -DMUMPS_LIBRARY_DMUMPS=/opt/homebrew/opt/brewsci-mumps/lib/libdmumps.dylib \
  -DMUMPS_LIBRARY_COMMON=/opt/homebrew/opt/brewsci-mumps/lib/libmumps_common.dylib \
  -DMUMPS_LIBRARY_PORD=/opt/homebrew/opt/brewsci-mumps/lib/libpord.dylib \
  -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"
```

3. **Build and install**:
```bash
make -j$(sysctl -n hw.ncpu)
sudo make install
```

4. **Fix library install names**:
```bash
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib
```

### Python Bindings Setup

1. **Fix Python module installation**:
```bash
sudo cp /usr/local/lib/python3.10/site-packages/akantu/akantu/__init__.py \
        /usr/local/lib/python3.10/site-packages/akantu/__init__.py
```

2. **Configure Python path** (add to `~/.zshrc`):
```bash
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"
```

3. **Remove mpi4py conflicts**:
```bash
conda remove mpi4py  # If present
```

4. **Verify installation**:
```bash
python -c "import akantu; print(akantu.__file__); print(len(dir(akantu)))"
```

Expected output: `/usr/local/lib/python3.10/site-packages/akantu/__init__.py` with approximately 221 attributes.

## Key Modifications from Original

### 1. macOS Build Compatibility (Required)

Two source files were modified to enable building on macOS Apple Silicon:

- **`cmake/Modules/FindMumps.cmake`**: Added macOS-specific workaround for MUMPS detection with MPI-enabled homebrew installation
- **`src/model/model.hh`**: Added conditional compilation guard for `CouplerSolidCohesiveContactOptions`

### 2. Phase Field Model Enhancements (Research Features)

Implemented tension-compression asymmetry in phase field fracture models:

- Added `phi_minus` and `phi_plus` internal fields for split energy contributions
- Implemented `computePhiMinusOnQuad()` in volumetric-deviatoric energy split
- Extended Python bindings for accessing asymmetric energy fields
- Modified quadratic phase field model to compute split energies

**Scientific motivation**: These enhancements enable modeling of materials with different behavior under tension and compression (e.g., concrete, ceramics, rocks), where cracks open under tension but do not interpenetrate under compression.

## Documentation

- **[BUILD_MACOS_README.md](BUILD_MACOS_README.md)**: Comprehensive build guide with troubleshooting
- **[CREATE_BACKUP.md](CREATE_BACKUP.md)**: Backup and restoration strategies
- **[ATTRIBUTION.md](ATTRIBUTION.md)**: Detailed attribution and modification documentation

## Known Issues and Solutions

### MUMPS Sequential Version Not Available
Homebrew's MUMPS is compiled with MPI support, not sequential. **Solution**: Use `AKANTU_PARALLEL=ON`.

### Scotch Linking Errors
Missing `_SCOTCH_errorPrint` symbol. **Solution**: Specify all three Scotch libraries in CMake configuration.

### Python Import Segfault
Conflict between mpi4py and Akantu's MPI initialization. **Solution**: Remove mpi4py from environment.

### Module Structure Issues
CMake installs `__init__.py` in nested directory. **Solution**: Copy to correct location (see Python Bindings Setup).

## Citation

If you use Akantu in your research, please cite the original work:

```bibtex
@article{akantu2024,
  title={Akantu: an HPC finite-element library for contact and dynamic fracture simulations},
  author={Richart, Nicolas and Molinari, Jean-Fran{\c{c}}ois and others},
  journal={The Journal of Open Source Software},
  year={2024},
  publisher={The Open Journal}
}
```

## License

This modified version maintains the original LGPLv3 license from Akantu.

```
Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)

Akantu is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.
```

Full license: http://www.gnu.org/licenses/

## Original Akantu Documentation

For general Akantu documentation, tutorials, and examples, please refer to:
- **Documentation**: https://akantu.readthedocs.io
- **Original Repository**: https://gitlab.com/akantu/akantu
- **Issues**: https://gitlab.com/akantu/akantu/-/issues

## Contributing

For issues specific to this macOS build, please open an issue on this repository. For general Akantu questions or contributions to the upstream project, please refer to the original GitLab repository.
