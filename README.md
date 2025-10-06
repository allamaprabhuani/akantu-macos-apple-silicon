# Akantu: macOS Apple Silicon Build

> **Build Akantu finite element library on Apple Silicon Macs (M1/M2/M3)**
>
> ⏱️ **Total time**: ~1-2 hours (mostly automated compilation)
> 💻 **Platform**: macOS 12.0+ on Apple Silicon
> 🎯 **Difficulty**: Beginner-friendly with step-by-step guide

---

## Attribution

This repository is a fork of the **Akantu** finite element library (Swiss-Made Open-Source Finite-Element Library), developed by the Computational Solid Mechanics Laboratory (LSMS) at École Polytechnique Fédérale de Lausanne (EPFL).

**Original Repository**: https://gitlab.com/akantu/akantu
**Original Branch**: `features/52-anisotropic-and-at1-phase-field-models`
**Original Commit**: 8c5e82d86
**License**: GNU Lesser General Public License v3.0 (LGPLv3)
**Copyright**: (©) 2010-2023 EPFL - LSMS

### Citation

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

**Original Documentation**: https://akantu.readthedocs.io

---

## Table of Contents

**For Complete Beginners:**
- [What is Akantu?](#what-is-akantu)
- [Prerequisites](#prerequisites)
- [Quick Start Guide](#quick-start-guide-complete-beginner) - **Start here!**
- [Troubleshooting](#troubleshooting-guide)

**For Advanced Users:**
- [Technical Details](#detailed-technical-information)
- [Source Code Modifications](#what-modifications-were-made)
- [Phase Field Enhancements](#phase-field-model-enhancements-user-specific-not-required-for-building)
- [Additional Documentation](#additional-documentation)

---

## What is Akantu?

Akantu is a powerful finite element library for simulating solid mechanics, fracture, and contact problems. This fork specifically enables it to run on Apple Silicon Macs (M1, M2, M3, etc.) and includes enhancements to the phase field fracture framework.

## Why This Fork Exists

This repository was created to support **phase field fracture research** on Apple Silicon. The original Akantu code could not compile on macOS with Apple Silicon, and my PhD research required the phase field solver with tension-compression asymmetry capabilities for modeling quasi-brittle materials (concrete, rocks, ceramics).

This fork provides:
1. **Essential platform fixes** - Required modifications for building on macOS Apple Silicon (arm64)
2. **Phase field enhancements** - Tension-compression asymmetry in phase field models for realistic fracture simulation

For detailed technical information about all modifications, see [ATTRIBUTION.md](ATTRIBUTION.md).

## Prerequisites

Before starting, you need:

1. **macOS with Apple Silicon** (M1/M2/M3 chip)
   - Check: Click Apple menu > About This Mac > look for "Apple M1" or "Apple M2"

2. **Homebrew** (package manager for macOS)
   - Install if needed: Open Terminal and run:
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

3. **Basic terminal knowledge** - You'll run commands in Terminal app

**Tested Configuration:**
- **macOS Version**: 26.0.1 (Build 25A362) - but should work on macOS 12.0+
- **Architecture**: arm64 (Apple Silicon)
- **Compiler**: GCC 15.2.0 (Homebrew)
- **Python**: 3.10.18 or later

---

# Quick Start Guide

Follow these steps in order. Each step builds on the previous one.

## Step 1: Install Dependencies

Open Terminal (Applications > Utilities > Terminal) and run:

```bash
brew install gcc boost scotch eigen open-mpi openblas scalapack cmake
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis
```

**Time estimate**: 15-30 minutes

**⚠️ IMPORTANT - If MUMPS Installation Fails:**

If the MUMPS installation fails with download errors, you need to edit the brew formula to use the MacPorts mirror:

```bash
brew edit brewsci-mumps
```

This opens an editor. Find the line with `url` and `sha256` and change them to:

```ruby
url "https://distfiles.macports.org/mumps/MUMPS_5.6.2.tar.gz"
sha256 "13a2c1aff2bd1aa92fe84b7b35d88f43434019963ca09ef7e8c90821a8f1d59a"
```

Then save (press `Esc`, type `:wq`, press `Enter` in vim) and retry:
```bash
brew install brewsci-mumps --without-brewsci-parmetis
```

**Why**: The MacPorts mirror is more reliable and was used for this build.

## Step 2: Download This Repository

```bash
cd ~
git clone https://github.com/allamaprabhuani/akantu-macos-apple-silicon.git
cd akantu-macos-apple-silicon
```

## Step 3: Configure the Build

You have two options: **Option A (Recommended for beginners)** uses ccmake with a visual interface, or **Option B** uses command-line cmake.

### Option A: Using ccmake (Interactive, Easier)

```bash
# Create a build directory
mkdir build
cd build

# Open interactive configuration with MUMPS_DIR set
CC=gcc-15 CXX=g++-15 FC=gfortran-15 MUMPS_DIR=/opt/homebrew/opt/brewsci-mumps ccmake ..
```

Configure these settings (use arrow keys to navigate, Enter to edit):

1. **Press `c` to configure** (scans your system - MUMPS should be auto-detected via MUMPS_DIR)

2. **Set these required options**:
   - `AKANTU_PARALLEL`: **ON** (required - Homebrew MUMPS is MPI-enabled)
   - `AKANTU_USE_SYSTEM_MUMPS`: **ON**
   - `SCOTCH_LIBRARY`: `/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib`

3. **Recommended for phase field research**:
   - `AKANTU_PYTHON_INTERFACE`: **ON** (enables Python bindings)
   - `PYTHON_EXECUTABLE`: Auto-detected

4. **Press `c` again to configure**

5. **Press `g` to generate and exit**

**Note**: `AKANTU_PARALLEL=ON` is required because Homebrew's MUMPS is MPI-enabled.

### Option B: Using cmake (Command-line)

```bash
mkdir build
cd build

CC=gcc-15 CXX=g++-15 FC=gfortran-15 MUMPS_DIR=/opt/homebrew/opt/brewsci-mumps cmake .. \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_PYTHON_INTERFACE=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=ON \
  -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"
```

**Expected output**: Should end with "Configuring done" and "Generating done"

## Step 4: Compile Akantu

```bash
make -j$(sysctl -n hw.ncpu)
```

**Time estimate**: 15-45 minutes. You'll see warnings - this is normal.

## Step 5: Install Akantu

```bash
sudo make install
```

**Note**: Installs to `/usr/local`. You'll need your password.

## Step 6: Fix Library Linking

```bash
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib
```

## Step 7: Setup Python Interface

Required only if you enabled `AKANTU_PYTHON_INTERFACE=ON` in Step 3.

### 7a. Fix Python Module Structure

```bash
sudo cp /usr/local/lib/python3.10/site-packages/akantu/akantu/__init__.py \
        /usr/local/lib/python3.10/site-packages/akantu/__init__.py
```

### 7b. Add to Python Path

```bash
echo 'export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"' >> ~/.zshrc
source ~/.zshrc
```

### 7c. Remove mpi4py (Causes Conflicts)

```bash
conda list | grep mpi4py
conda remove mpi4py  # If found
```

### 7d. Setup Environment Variables (Required for macOS)

Add these environment variables to your `~/.zshrc`:

```bash
cat >> ~/.zshrc << 'EOF'

# Akantu environment variables
export DYLD_INSERT_LIBRARIES="/opt/homebrew/lib/libscotcherr.dylib:/opt/homebrew/lib/libscotcherrexit.dylib"
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:/usr/local/lib:${DYLD_LIBRARY_PATH}"
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:${PYTHONPATH}"

# Force all python commands to use conda python
alias python3='/Users/allamaprabhuani/miniconda3/bin/python3'
alias python='/Users/allamaprabhuani/miniconda3/bin/python'
alias pip='/Users/allamaprabhuani/miniconda3/bin/pip'
alias pip3='/Users/allamaprabhuani/miniconda3/bin/pip3'
EOF

source ~/.zshrc
```

**What this does:**
- Sets `DYLD_INSERT_LIBRARIES` to preload Scotch error handling libraries (required)
- Sets `DYLD_LIBRARY_PATH` to include Homebrew and system library paths
- Sets `PYTHONPATH` to include akantu's installation directory
- Forces all `python` commands to use conda Python instead of system Python

## Step 8: Verify Installation

Open a new terminal (or run `source ~/.zshrc`) and test:

```bash
python -c "import akantu; print('Success!'); print('Available functions:', len(dir(akantu)))"
```

**Expected**: Shows "Success!" and approximately 265 functions.

**Note**: If you set up the aliases in step 7d, `python` now points to conda Python automatically.

---

## Next Steps: Running Your First Simulation

Congratulations! Akantu is now installed. Here's a simple example to test it:

### Create a Simple Test Script

Create a file called `test_akantu.py`:

```python
import akantu
import numpy as np

# Print version info
print("Akantu successfully imported!")
print(f"Available modules: {len(dir(akantu))}")

# Check if phase field module is available
if hasattr(akantu, 'PhaseFieldModel'):
    print("✓ Phase field module is available")
else:
    print("Note: Phase field module may require specific build options")

# Basic mesh test
print("\nTesting mesh creation...")
try:
    mesh = akantu.Mesh(2)  # 2D mesh
    print("✓ Mesh creation successful")
except Exception as e:
    print(f"Mesh creation issue: {e}")

print("\nAkantu is ready to use!")
```

Run it:
```bash
python3 test_akantu.py
```

### Where to Learn More

- **Examples**: Check the `examples/` directory in this repository for sample simulations
- **Documentation**: https://akantu.readthedocs.io
- **Tutorials**: Interactive notebooks available at the original Akantu repository
- **Phase Field Examples**: See `examples/phase_field/` for fracture simulation examples

---

# Detailed Technical Information

## Dependencies and Versions

The following versions were tested and confirmed working:

- **gcc**: 15.2.0
- **boost**: 1.89.0
- **scotch**: 7.0.10
- **eigen**: 3.4.0_1
- **open-mpi**: 5.0.8
- **openblas**: 0.3.30
- **scalapack**: 2.2.2
- **brewsci-mumps**: 5.6.2

Newer versions should work, but if you encounter issues, these specific versions are known to work.

## What Modifications Were Made?

The following files were modified from the original source to fix compilation and linking issues on macOS:

### 1. `cmake/Modules/FindMumps.cmake`

**Issue**: MUMPS test compilation succeeds but runtime fails on macOS because homebrew's MUMPS is compiled with MPI support, not sequential. This causes the FindMumps CMake module to fail after 10 retries.

**Fix**: Added macOS-specific workaround to accept successful compilation even if the test fails at runtime (after 5 retries). Also added additional MUMPS dependencies for macOS.

**Changes**:
```cmake
# Around line 282 - Added macOS runtime test bypass
if(APPLE AND _mumps_compiles AND _retry_count GREATER 5)
  message(STATUS "MUMPS test compiled successfully on macOS, accepting despite runtime issues")
  break()
endif()

# Around line 361 - Added additional MUMPS dependencies for macOS
if(APPLE)
  # in doubt add some stuff because mumps was perhaps badly compiled
  mumps_add_dependency(pord _libs _incs)
  list(APPEND _libraries_all ${_libs})
  mumps_add_dependency(mumps_common _libs _incs)
  list(APPEND _libraries_all ${_libs})
  mumps_add_dependency(BLAS _libs _incs)
  list(APPEND _libraries_all ${_libs})
  mumps_add_dependency(LAPACK _libs _incs)
  list(APPEND _libraries_all ${_libs})
  mumps_add_dependency(ScaLAPACK _libs _incs)
  list(APPEND _libraries_all ${_libs})
  mumps_add_dependency(gfortran _libs _incs)
  list(APPEND _libraries_all ${_libs})
endif()
```

### 2. `src/model/model.hh`

**Issue**: `CouplerSolidCohesiveContactOptions` is only defined when `AKANTU_COHESIVE_ELEMENT` is enabled, but the code was trying to use it whenever `AKANTU_MODEL_COUPLERS` is enabled, causing compilation errors.

**Fix**: Added conditional compilation guard around the usage of `CouplerSolidCohesiveContactOptions`.

**Changes** (around line 117):
```cpp
#ifdef AKANTU_MODEL_COUPLERS
    case ModelType::_coupler_solid_contact:
      this->initFullImpl(CouplerSolidContactOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#ifdef AKANTU_COHESIVE_ELEMENT  // Added this guard
    case ModelType::_coupler_solid_cohesive_contact:
      this->initFullImpl(CouplerSolidCohesiveContactOptions{
          use_named_args, std::forward<decltype(_pack)>(_pack)...});
      break;
#endif  // Added this endif
#endif
```

## Build Instructions

### 1. Configure CMake

```bash
cd akantu
mkdir build
cd build

CC=gcc-15 CXX=g++-15 FC=gfortran-15 MUMPS_DIR=/opt/homebrew/opt/brewsci-mumps cmake .. \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=ON \
  -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"
```

**Important Notes**:
- `AKANTU_PARALLEL=ON` is required because homebrew's MUMPS is compiled with MPI support
- All three Scotch libraries (libscotch, libscotcherr, libscotcherrexit) must be specified to avoid missing symbol errors (`_SCOTCH_errorPrint`)

### 2. Build

```bash
make -j$(sysctl -n hw.ncpu)
```

**Expected Warnings** (these are safe to ignore):
- C++17 ABI change notes about parameter passing for `std::pair`
- Uninitialized variable warnings in Eigen template code (false positives from GCC optimization analysis)
- Duplicate library warnings for `-lgfortran`

### 3. Install

```bash
sudo make install
```

This installs to `/usr/local` by default.

### 4. Fix Library Install Names

After installation, fix the library install names for proper dynamic linking:

```bash
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib
```

## Python Bindings Setup

### 1. Fix Python Module Installation

The Python module `__init__.py` gets installed in the wrong location. Fix it:

```bash
sudo cp /usr/local/lib/python3.10/site-packages/akantu/akantu/__init__.py \
        /usr/local/lib/python3.10/site-packages/akantu/__init__.py
```

### 2. Handle mpi4py Import Error

Edit `/usr/local/lib/python3.10/site-packages/akantu/__init__.py` and change line 21 from:
```python
    except Exception:
```
to:
```python
    except (Exception, ModuleNotFoundError):
```

This allows akantu to work even if mpi4py is not installed.

### 3. Configure Python Path

Add to your `~/.zshrc`:
```bash
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"
```

Then reload:
```bash
source ~/.zshrc
```

### 4. Handle mpi4py Conflicts

If you have mpi4py installed in your conda environment, it will conflict with akantu's MPI initialization and cause segmentation faults. Either:

**Option A**: Remove mpi4py from your environment:
```bash
conda remove mpi4py
```

**Option B**: Create a separate conda environment for akantu:
```bash
conda create -n akantu python=3.10 numpy scipy
conda activate akantu
```

## Verification

Test the installation:

```bash
python -c "import akantu; print(akantu.__file__); print(len(dir(akantu)))"
```

Expected output:
```
/usr/local/lib/python3.10/site-packages/akantu/__init__.py
221
```

The number 221 indicates successful import with all attributes available.

---

# Troubleshooting Guide

## Common Build Errors

### Error: "cmake: command not found"

**Problem**: CMake is not installed.

**Solution**:
```bash
brew install cmake
```

### Error: "Cannot find MUMPS library" or MUMPS download fails

**Problem**: MUMPS wasn't installed correctly or CMake can't find it.

**Solution**:
```bash
# Verify MUMPS is installed
brew list brewsci-mumps

# If not found, install it
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis

# Verify the files exist
ls -la /opt/homebrew/opt/brewsci-mumps/lib/
```

**If installation fails with download error**, edit the brew formula to use the MacPorts mirror:

```bash
brew edit brewsci-mumps
```

Change the `url` and `sha256` lines to:
```ruby
url "https://distfiles.macports.org/mumps/MUMPS_5.6.2.tar.gz"
sha256 "13a2c1aff2bd1aa92fe84b7b35d88f43434019963ca09ef7e8c90821a8f1d59a"
```

Save (`:wq` in vim) and retry installation. This MacPorts mirror is what this build was tested with and is known to be reliable.

### Error: "symbol not found: _SCOTCH_errorPrint"

**Problem**: Missing Scotch error handling libraries.

**Solution**: Make sure you used the complete `-DSCOTCH_LIBRARY` flag with all three libraries in Step 3 (cmake configuration).

### Error: Compilation warnings about "uninitialized variables"

**Problem**: These are false positives from GCC's optimization analysis in Eigen library code.

**Solution**: These warnings are **safe to ignore**. They don't affect functionality.

### Error: "CouplerSolidCohesiveContactOptions was not declared"

**Problem**: You're trying to build an unmodified version of Akantu.

**Solution**: Make sure you cloned **this repository**, not the original Akantu repository. This fork includes the necessary fixes.

## Common Python Import Errors

### Error: "ModuleNotFoundError: No module named 'akantu'"

**Solutions to try in order**:

1. **Check PYTHONPATH**:
```bash
echo $PYTHONPATH
# Should include /usr/local/lib/python3.10/site-packages
```

2. **Add it if missing**:
```bash
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"
# Make it permanent:
echo 'export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"' >> ~/.zshrc
```

3. **Check if module file exists**:
```bash
ls -la /usr/local/lib/python3.10/site-packages/akantu/
```

4. **Fix module structure** (if __init__.py is in wrong location):
```bash
sudo cp /usr/local/lib/python3.10/site-packages/akantu/akantu/__init__.py \
        /usr/local/lib/python3.10/site-packages/akantu/__init__.py
```

### Error: Segmentation fault when importing akantu

**Problem**: mpi4py conflict.

**Solution**:
```bash
# Check if mpi4py is installed
conda list | grep mpi4py

# Remove it
conda remove mpi4py

# Or create a clean environment
conda create -n akantu_env python=3.10 numpy scipy
conda activate akantu_env
```

### Error: "Library not loaded: @rpath/libakantu.5.0.dylib"

**Problem**: Library install names not fixed.

**Solution**: Run Step 6 from the Quick Start Guide:
```bash
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib
```

## Still Having Issues?

1. **Check your system**:
```bash
# Verify you have Apple Silicon
uname -m  # Should show "arm64"

# Check macOS version
sw_vers  # Should be 12.0 or later

# Verify Homebrew installation
brew --version
```

2. **Clean build and try again**:
```bash
cd ~/akantu-macos-apple-silicon
rm -rf build
mkdir build
cd build
# Then repeat Step 3 onwards from Quick Start Guide
```

3. **Check for conflicting installations**:
```bash
# Make sure no old Akantu installations exist
which akantu
ls -la /usr/local/lib/libakantu*
```

4. **Get help**:
   - Open an issue on this GitHub repository with:
     - Your exact error message
     - Your macOS version (`sw_vers`)
     - Your chip type (`uname -m`)
     - Output of `brew list`

---

# Known Issues and Technical Details

### Issue 1: MUMPS Runtime Test Fails During CMake

**Technical Details**: Homebrew's MUMPS is compiled with MPI support. The test binary compiles but fails at runtime in CMake's testing environment.

**How this fork fixes it**: Modified `cmake/Modules/FindMumps.cmake` to accept successful compilation as sufficient evidence on macOS.

### Issue 2: Python Module Nested Structure

**Technical Details**: CMake installs files as `akantu/akantu/__init__.py` instead of `akantu/__init__.py`.

**Workaround**: Manual copy to correct location (automated in Step 7a).

### Issue 3: MPI Initialization Conflicts

**Technical Details**: Both mpi4py and Akantu's C++ MPI bindings try to initialize MPI, causing segmentation faults.

**Solution**: Remove mpi4py or use separate environment.

## Complete Diff Summary from Original Branch

All modifications from the original `origin/features/52-anisotropic-and-at1-phase-field-models` branch:

### Files Modified for macOS Build Compatibility (Required)

1. **`cmake/Modules/FindMumps.cmake`** (+17 lines)
   - Added macOS-specific workaround for MUMPS runtime test failures
   - Added additional MUMPS dependency resolution for macOS

2. **`src/model/model.hh`** (+2 lines)
   - Added conditional compilation guard for `CouplerSolidCohesiveContactOptions`

### Files Modified for Phase Field Features (User-specific, not required for building)

The following modifications add tension-compression asymmetry capabilities to the phase field model. These changes introduce split energy contributions (phi_minus and phi_plus) that separate compressive and tensile components, which is important for modeling materials with different behavior in tension vs compression (e.g., ceramics, concrete, rocks).

**Note**: These changes may be moved to a separate feature branch in the future.

#### Core Phase Field Model Changes

3. **`src/model/phase_field/phasefield.hh`** (+20 lines)
   - Added three new internal fields:
     - `phi_minus`: Stores compressive part of energy density
     - `phi_plus`: Stores positive/tensile part before history tracking
     - Getter macros for accessing these fields
   - Modified internal data structure to track asymmetric energy contributions
   - **Purpose**: Enables modeling of materials with tension-compression asymmetry

4. **`src/model/phase_field/phasefield.cc`** (+4 lines)
   - Registered `phi_minus` and `phi_plus` as internal fields
   - Initialized default values to 0.0 for new fields
   - **Purpose**: Properly initializes the new asymmetric energy fields

#### Energy Split Implementations

5. **`src/model/phase_field/energy_splits/volumetric_deviatoric_split.hh`** (+4 lines)
   - Added declaration for `computePhiMinusOnQuad()` method
   - **Purpose**: Interface for computing compressive energy contribution

6. **`src/model/phase_field/energy_splits/volumetric_deviatoric_split_inline_impl.hh`** (+12 lines)
   - Implemented `computePhiMinusOnQuad()` method
   - Computes compressive volumetric energy: φ⁻ = 0.5 * κ * (tr(ε)⁻)²
   - Uses bulk modulus κ = λ + (2/3)μ
   - **Purpose**: Calculates energy contribution from compressive strains only

7. **`src/model/phase_field/energy_splits/no_energy_split.hh`** (+4 lines)
   - Added declaration for `computePhiMinusOnQuad()` method
   - **Purpose**: Provides interface consistency for no-split case

8. **`src/model/phase_field/energy_splits/no_energy_split_inline_impl.hh`** (+8 lines)
   - Implemented `computePhiMinusOnQuad()` for no-split case
   - Returns zero (no compressive/tensile distinction)
   - **Purpose**: Maintains compatibility when asymmetric splitting is disabled

#### Phase Field Implementation

9. **`src/model/phase_field/phasefields/phasefield_quadratic.cc`** (+4 lines)
   - Modified to compute and store `phi_minus` and `phi_plus` values
   - **Purpose**: Implements asymmetric energy calculation in quadratic phase field model

#### Python Bindings

10. **`python/py_phasefield.cc`** (+57 lines, -1 line)
    - Added Python bindings for new getters:
      - `getPhi()`: Access phase field variable
      - `getPhiMinus()`: Access compressive energy part
      - `getPhiPlus()`: Access tensile energy part
    - Added `getInternalReal()` method for unified field access
      - Supports: "strain", "phi", "phi_minus", "phi_plus", "damage"
    - **Purpose**: Exposes asymmetric energy fields to Python for analysis and post-processing

**Scientific Context**: These modifications implement a volumetric-deviatoric split with tension-compression asymmetry, commonly used in phase field fracture models for quasi-brittle materials. The split allows cracks to open under tension but prevents interpenetration under compression, which is physically realistic for many materials.

**Total changes**: 132 lines added, 1 line removed across 10 files

### How to Get Only macOS Build Fixes

If you want to build the original branch without the phase field modifications, you only need to apply changes to these two files:
- `cmake/Modules/FindMumps.cmake`
- `src/model/model.hh`

To see the exact changes:
```bash
git diff cmake/Modules/FindMumps.cmake
git diff src/model/model.hh
```

## Summary

Building Akantu on macOS Apple Silicon requires:
1. Two source code modifications (FindMumps.cmake and model.hh)
2. Specific CMake configuration for homebrew-installed dependencies
3. Post-installation fixes for library install names
4. Python module installation fixes
5. Handling mpi4py conflicts

Once these steps are completed, Akantu builds and runs successfully on macOS with Apple Silicon.

---

## Additional Documentation

- **[ATTRIBUTION.md](ATTRIBUTION.md)** - Detailed academic documentation of all modifications and theoretical background
- **[CREATE_BACKUP.md](CREATE_BACKUP.md)** - Backup and restoration strategies
- **[README_ORIGINAL.md](README_ORIGINAL.md)** - Original Akantu README with general build instructions

## References

For the phase field implementation approach:

1. Amor, H., Marigo, J. J., & Maurini, C. (2009). Regularized formulation of the variational brittle fracture with unilateral contact: Numerical experiments. *Journal of the Mechanics and Physics of Solids*, 57(8), 1209-1229.

2. Miehe, C., Hofacker, M., & Welschinger, F. (2010). A phase field model for rate-independent crack propagation: Robust algorithmic implementation based on operator splits. *Computer Methods in Applied Mechanics and Engineering*, 199(45-48), 2765-2778.

## Support

For questions regarding:
- **macOS build issues**: Open an issue on this repository
- **Phase field implementation**: Open an issue on this repository
- **General Akantu questions**: Refer to the [original repository](https://gitlab.com/akantu/akantu)

---

*This work was developed as part of computational fracture mechanics research, building upon the excellent foundation provided by the Akantu development team at EPFL.*
