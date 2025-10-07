# Akantu Wheel Technical Documentation

## What is the Wheel?

Pre-compiled Python package for Akantu v5.0.7 on macOS Apple Silicon.

**File**: `akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl` (38.6 MB)
**Download**: [GitHub Releases](https://github.com/allamaprabhuani/akantu-macos-apple-silicon/releases/tag/v5.0.7-macos-wheel)

## Contents

- Python extension: `py11_akantu.cpython-310-darwin.so`
- Akantu libraries: `libakantu.5.0.dylib`, `libiohelper.dylib`
- 26 bundled dependencies (MPI, Scotch, MUMPS, GCC runtime, math libraries)

All libraries use `@loader_path` references for portability.

## Platform Compatibility

**Works on:**
- macOS 11.0+ (Big Sur or later)
- Apple Silicon (M1/M2/M3/M4)
- Python 3.10 only

**Does not work on:**
- Intel Macs (x86_64)
- Other Python versions (3.11, 3.12, etc.)
- Linux or Windows

## Python Version Specificity

The wheel is compiled for Python 3.10 specifically. The C extension (`cpython-310`) is ABI-incompatible with other Python versions.

**Check your version:**
```bash
python3 --version  # Must show 3.10.x
```

**For other versions:**
```bash
python3 packaging/create_standalone_wheel.py
```

## Wheel Naming Convention

`akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl`

- **akantu** - Package name
- **5.0.7** - Version
- **cp310** - CPython 3.10
- **none** - No specific ABI
- **macosx_26_0_arm64** - macOS 26.0+, Apple Silicon

## Dependencies

### Bundled in Wheel
- MPI: Open-MPI libraries
- Solvers: MUMPS, ScaLAPACK, OpenBLAS
- Partitioning: Scotch libraries
- Runtime: GCC (gfortran, libgomp, libquadmath)
- System: hwloc, libevent, lzma

### Required on System
Homebrew packages needed for symbol resolution:

```bash
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis
```

**Why?** Scotch uses weak symbols across libraries that must be resolved at runtime via system libraries.

## Build Configuration

```cmake
AKANTU_PARALLEL=ON
AKANTU_PYTHON_INTERFACE=ON
AKANTU_USE_SYSTEM_MUMPS=ON
SCOTCH_LIBRARY=libscotch+libscotcherr+libscotcherrexit
```

Compiler: GCC 15.2.0

## Creating Custom Wheels

### For Different Python Version

```bash
python3.11 packaging/create_standalone_wheel.py
# Creates: akantu-5.0.7-cp311-none-macosx_26_0_arm64.whl
```

### Scripts

- **`create_standalone_wheel.py`** - Bundles all dependencies (recommended)
- **`create_akantu_wheel_mac.py`** - Minimal wheel requiring system libs
- **`test_wheel_install.sh`** - Tests installation

### Build Process

1. Recursively collects all `.dylib` dependencies
2. Copies to `.dylibs/` directory in package
3. Fixes library paths using `install_name_tool`
4. Modifies `__init__.py` to preload Scotch error libraries
5. Creates wheel with correct platform tags

## Known Limitations

1. **Requires Homebrew**: Despite bundling libraries, system Scotch needed for symbol resolution
2. **Python 3.10 only**: Binary incompatibility across Python versions
3. **Apple Silicon only**: Architecture-specific compilation
4. **Not truly standalone**: Dependency on system libraries for weak symbols

## Verification

**SHA256 Checksum:**
```bash
shasum -a 256 akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl
# Should match: 3e2d0d9f25b041883ad38718112158c4f50bc0bd41828a354b3a7ca576ae7d5e
```

See `packaging/CHECKSUM.txt`

## Troubleshooting

### Symbol not found errors
Install Homebrew dependencies, particularly Scotch.

### Import fails completely
- Check Python version matches wheel (3.10)
- Verify dependencies installed: `brew list scotch`

### Library path errors
Check bundled libraries exist:
```bash
unzip -l akantu-5.0.7-*.whl | grep dylibs
```

### Works in one environment but not another
Virtual environments must have same Python 3.10 version and Homebrew dependencies.

## Performance

| Aspect | Source Build | Wheel Install |
|--------|-------------|---------------|
| Time | 1-2 hours | 10 seconds |
| Disk during build | ~5 GB | ~40 MB |
| Result size | ~500 MB | ~40 MB |
| Customization | Full | None |
| Portability | One machine | Similar configs |

## License Compliance

Wheel bundles libraries under various licenses:
- Akantu: LGPLv3
- GCC runtime: GPLv3 with runtime exception
- Scotch: CeCILL-C (LGPL-compatible)
- MUMPS: CeCILL-C
- OpenBLAS: BSD
- Open-MPI: BSD

Distribution requires source code availability (LGPLv3).

## For Developers

### Inspect wheel contents
```bash
unzip -l akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl
```

### Check library dependencies
```bash
otool -L /path/to/site-packages/akantu/py11_akantu.cpython-310-darwin.so
```

### Debug symbol resolution
```bash
nm -g /path/to/site-packages/akantu/.dylibs/libscotch.7.0.dylib | grep SCOTCH_errorPrint
```

## Future Improvements

1. Static linking of Scotch to eliminate system dependency
2. Multi-Python version support via delocate tool
3. Automated CI/CD for wheel building
4. Universal binary for arm64 and x86_64

## References

- Wheel specification: [PEP 427](https://www.python.org/dev/peps/pep-0427/)
- Platform tags: [PEP 425](https://www.python.org/dev/peps/pep-0425/)
- Official docs: https://akantu.readthedocs.io
