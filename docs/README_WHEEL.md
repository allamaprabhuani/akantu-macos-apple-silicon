# Akantu Wheel Package for macOS Apple Silicon

This directory contains scripts and wheels for packaging your built Akantu installation as Python wheels for easy installation.

## Overview

You've successfully built Akantu from source on your Mac (which takes 1-2 hours). Now you can package it as a wheel to:
- Install quickly on this Mac without rebuilding
- Transfer to other similar Macs
- Use in virtual environments
- Manage with pip like any Python package

## Files

- **`create_standalone_wheel.py`** - Creates a wheel with all dynamic libraries bundled (~39 MB)
- **`create_akantu_wheel_mac.py`** - Creates a smaller wheel that requires Homebrew dependencies (~16 MB)
- **`test_wheel_install.sh`** - Tests wheel installation in a clean environment
- **`dist/akantu-5.0.7-py3-none-any.whl`** - The generated wheel file

## Quick Start

### Option 1: Use the Wheel on This Mac

The wheel has already been created. Just install it:

```bash
pip install --user dist/akantu-5.0.7-py3-none-any.whl
```

Then verify:

```bash
python3 -c "import akantu; print('Success!')"
```

**Note**: This requires that you have the same Homebrew dependencies installed that were used during the build:

```bash
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis
```

### Option 2: Transfer to Another Mac

The wheel will work on other Macs with:
- Apple Silicon (M1/M2/M3)
- macOS 11.0 or later
- Python 3.10
- Same Homebrew dependencies installed

Transfer the wheel:

```bash
scp dist/akantu-5.0.7-py3-none-any.whl user@othermac:
```

On the other Mac:

```bash
# Install Homebrew dependencies first
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis

# Install the wheel
pip install --user akantu-5.0.7-py3-none-any.whl
```

## Rebuild the Wheel

If you update your Akantu build and want to recreate the wheel:

```bash
# Make sure Akantu is built
cd akantu/build
make -j$(sysctl -n hw.ncpu)

# Return to parent directory and create wheel
cd ../..
python3 create_standalone_wheel.py
```

This will create a new wheel in `dist/`.

## Technical Details

### What's in the Wheel?

The standalone wheel contains:

1. **Python Extension**: `py11_akantu.cpython-310-darwin.so` (9.5 MB)
2. **Akantu Libraries**:
   - `libakantu.5.0.dylib` (~114 MB)
   - `libiohelper.dylib` (~519 KB)
3. **Dependencies** (26 libraries total, ~39 MB):
   - MPI libraries (Open-MPI)
   - Scotch libraries (graph partitioning)
   - MUMPS libraries (sparse solver)
   - GCC runtime libraries (gfortran, gomp, etc.)
   - Math libraries (OpenBLAS, ScaLAPACK)
   - Utility libraries (hwloc, libevent, etc.)

### Library Path Handling

All libraries in the wheel use `@loader_path/.dylibs/` references, so they can find each other regardless of where the wheel is installed.

The creation script:
1. Recursively collects all `.dylib` dependencies
2. Copies them to `.dylibs/` directory
3. Fixes all library paths using `install_name_tool`
4. Packages everything into a single wheel file

### Limitations

**Important**: The wheel has Scotch library symbol resolution issues that require the Homebrew scotch libraries to be present. This is because:

1. Scotch uses weak symbols (`_SCOTCH_errorPrint`) across multiple libraries
2. These symbols must be resolved at runtime
3. macOS's dyld requires all dependencies to be explicitly linked or in the library path

**Current Status**: The wheel bundles all libraries but still requires Homebrew to be installed on the target system for symbol resolution to work properly.

## Alternative Approaches Tried

1. **Preloading with ctypes.CDLL**: Doesn't add symbols to global namespace early enough
2. **Using @rpath**: Requires complex RPATH management
3. **install_name_tool changes**: Successfully updates paths but doesn't resolve weak symbols
4. **Bundling all dependencies**: Works for most libraries except Scotch

## Workaround for Full Standalone

If you need a truly standalone version that works without Homebrew, consider:

1. **Static linking**: Rebuild Akantu with static libraries (requires rebuilding from source)
2. **Conda package**: Use conda-build to create a fully self-contained package
3. **Docker**: Package everything in a container

## Usage After Installation

```python
import akantu
import numpy as np

# Create a mesh
mesh = akantu.Mesh(spatial_dimension=2)
mesh.read("your_mesh.msh")

# Create a model
model = akantu.SolidMechanicsModel(mesh)
model.initFull()

# Your simulation code here...
```

For examples, see: https://akantu.readthedocs.io

## Troubleshooting

### ImportError: symbol not found '_SCOTCH_errorPrint'

This means the Scotch libraries aren't being found. Solutions:

1. Install Homebrew scotch: `brew install scotch`
2. Use your original build installation instead of the wheel
3. Set `DYLD_LIBRARY_PATH` (not recommended for security reasons)

### Wrong Python version

The wheel is built for Python 3.10. If you're using a different version:

```bash
# Check your version
python3 --version

# If different, rebuild the wheel with your Python version
python3 create_standalone_wheel.py
```

### Library loading errors

```bash
# Check what libraries are being loaded
otool -L /path/to/site-packages/akantu/py11_akantu.cpython-310-darwin.so

# Check if libraries exist
ls -lh /path/to/site-packages/akantu/.dylibs/
```

## Benefits vs. Building from Source

| Aspect | Building from Source | Using Wheel |
|--------|---------------------|-------------|
| Time | 1-2 hours | 10 seconds |
| Customization | Full control | Pre-built configuration |
| Dependencies | Must install | Must install (same) |
| Portability | One machine | Multiple machines |
| Updates | Rebuild everything | Just reinstall wheel |
| Virtual envs | System-wide only | Works in venvs |

## Future Improvements

To make the wheel truly standalone:

1. Statically link Scotch during Akantu build
2. Use `delocate` tool (like `auditwheel` for macOS)
3. Create a conda package instead of a wheel
4. Build universal binaries for both arm64 and x86_64

## License

Akantu is licensed under LGPLv3. When distributing the wheel, ensure compliance with:
- Akantu's LGPL license
- Licenses of bundled dependencies (GCC runtime, Scotch, MUMPS, etc.)

## Getting Help

- **Akantu Documentation**: https://akantu.readthedocs.io
- **Akantu Issues**: https://gitlab.com/akantu/akantu/-/issues
- **macOS Build Fork**: https://github.com/allamaprabhuani/akantu-macos-apple-silicon

## Summary

✅ **What Works**: The wheel successfully packages Akantu with all dependencies
✅ **Easy to Use**: Simple `pip install` command
✅ **Portable**: Transfer to similar Macs
⚠️ **Caveat**: Still requires Homebrew scotch libraries on target system

For a quick installation on this Mac or similar Macs with Homebrew, this wheel saves significant build time!
