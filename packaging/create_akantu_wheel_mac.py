#!/usr/bin/env python3
"""
Create a wheel from an already-built Akantu installation (Mac version)

This script packages the pre-built Akantu binaries into a wheel file
for easy installation via pip.

Usage:
    python3 create_akantu_wheel_mac.py

After running, install with:
    pip install dist/akantu-*.whl
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
import re

# Configuration
AKANTU_ROOT = Path(__file__).parent / "akantu"
BUILD_DIR = AKANTU_ROOT / "build"
PYTHON_BUILD_DIR = BUILD_DIR / "python" / "akantu"
WHEEL_BUILD_DIR = Path(__file__).parent / "wheel_build"
DIST_DIR = Path(__file__).parent / "dist"

def get_akantu_version():
    """Extract version from CMakeLists.txt or use a default"""
    cmake_file = AKANTU_ROOT / "CMakeLists.txt"
    if cmake_file.exists():
        content = cmake_file.read_text()
        # Look for version definition
        match = re.search(r'project\([^)]*VERSION\s+([0-9.]+)', content, re.IGNORECASE)
        if match:
            return match.group(1)

    # Fallback: check built files
    version_file = BUILD_DIR / "VERSION"
    if version_file.exists():
        return version_file.read_text().strip()

    return "5.0.7"  # Default version

def get_python_version():
    """Get Python version tag for wheel"""
    return f"cp{sys.version_info.major}{sys.version_info.minor}"

def get_platform_tag():
    """Get platform tag for wheel"""
    # For macOS: macosx_<major>_<minor>_<arch>
    import platform
    mac_ver = platform.mac_ver()[0].split('.')
    arch = platform.machine()  # 'arm64' or 'x86_64'

    # Convert to wheel format
    if arch == 'arm64':
        return f"macosx_{mac_ver[0]}_{mac_ver[1]}_arm64"
    else:
        return f"macosx_{mac_ver[0]}_{mac_ver[1]}_x86_64"

def create_wheel_structure():
    """Create the wheel directory structure"""
    print("=" * 60)
    print("Creating wheel structure")
    print("=" * 60)

    # Clean and create directories
    if WHEEL_BUILD_DIR.exists():
        shutil.rmtree(WHEEL_BUILD_DIR)
    WHEEL_BUILD_DIR.mkdir(parents=True)

    if DIST_DIR.exists():
        shutil.rmtree(DIST_DIR)
    DIST_DIR.mkdir(parents=True)

    # Create package directory
    package_dir = WHEEL_BUILD_DIR / "akantu"
    package_dir.mkdir()

    # Copy Python module files
    print(f"Copying from: {PYTHON_BUILD_DIR}")

    # Copy __init__.py
    init_src = PYTHON_BUILD_DIR / "__init__.py"
    init_dst = package_dir / "__init__.py"
    if init_src.exists():
        shutil.copy2(init_src, init_dst)
        print(f"  ✓ Copied __init__.py")
    else:
        print(f"  ✗ Warning: {init_src} not found")

    # Copy the compiled Python extension
    so_files = list(PYTHON_BUILD_DIR.glob("*.so"))
    if so_files:
        for so_file in so_files:
            shutil.copy2(so_file, package_dir / so_file.name)
            print(f"  ✓ Copied {so_file.name}")
    else:
        print(f"  ✗ Error: No .so files found in {PYTHON_BUILD_DIR}")
        return False

    # Copy dynamic libraries that the module depends on
    lib_dir = BUILD_DIR / "src"
    dylib_files = list(lib_dir.glob("*.dylib"))
    if dylib_files:
        libs_dir = package_dir / ".dylibs"
        libs_dir.mkdir()
        for dylib in dylib_files:
            # Only copy the main symlink and actual file
            if not dylib.is_symlink():
                shutil.copy2(dylib, libs_dir / dylib.name)
                print(f"  ✓ Copied {dylib.name}")

    # Copy iohelper library if present
    iohelper_lib = BUILD_DIR / "third-party" / "iohelper" / "src" / "libiohelper.dylib"
    if iohelper_lib.exists():
        if not (package_dir / ".dylibs").exists():
            (package_dir / ".dylibs").mkdir()
        shutil.copy2(iohelper_lib, package_dir / ".dylibs" / "libiohelper.dylib")
        print(f"  ✓ Copied libiohelper.dylib")

    return True

def create_setup_py():
    """Create a simple setup.py for wheel building"""
    setup_content = '''#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import os

# Get all package data (shared libraries)
package_data = []
akantu_dir = os.path.join(os.path.dirname(__file__), "akantu")
if os.path.exists(os.path.join(akantu_dir, ".dylibs")):
    package_data.append(".dylibs/*")
package_data.append("*.so")

setup(
    name="akantu",
    version="{version}",
    url="https://akantu.ch",
    author="Nicolas Richart",
    author_email="nicolas.richart@epfl.ch",
    description="Akantu: Swiss-Made Open-Source Finite-Element Library (Pre-built for Mac)",
    long_description="Pre-compiled Akantu library for macOS. Install with: pip install akantu-{version}-*.whl",
    long_description_content_type="text/plain",
    platforms="macOS",
    license="L-GPLv3",
    project_urls={{
        "Bug Tracker": "https://gitlab.com/akantu/akantu/-/issues",
        "Source": "https://gitlab.com/akantu/akantu",
    }},
    install_requires=["numpy", "scipy"],
    packages=find_packages(),
    package_data={{"akantu": package_data}},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Topic :: Education",
        "Topic :: Scientific/Engineering",
    ],
)
'''

    version = get_akantu_version()
    setup_path = WHEEL_BUILD_DIR / "setup.py"
    setup_path.write_text(setup_content.format(version=version))
    print(f"\n✓ Created setup.py with version {version}")
    return version

def build_wheel():
    """Build the wheel using setuptools"""
    print("\n" + "=" * 60)
    print("Building wheel")
    print("=" * 60)

    os.chdir(WHEEL_BUILD_DIR)

    # Build the wheel
    result = subprocess.run(
        [sys.executable, "setup.py", "bdist_wheel"],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("Error building wheel:")
        print(result.stderr)
        return False

    print(result.stdout)

    # Move wheel to dist directory
    wheel_dist = WHEEL_BUILD_DIR / "dist"
    if wheel_dist.exists():
        for wheel_file in wheel_dist.glob("*.whl"):
            dest = DIST_DIR / wheel_file.name
            shutil.copy2(wheel_file, dest)
            print(f"\n✓ Wheel created: {dest}")
            print(f"  Size: {dest.stat().st_size / (1024*1024):.1f} MB")

    return True

def create_installation_guide(version):
    """Create installation instructions"""
    import platform

    guide = f'''
================================================================================
Akantu Wheel Installation Guide (macOS)
================================================================================

WHEEL INFORMATION:
  Version: {version}
  Platform: {get_platform_tag()}
  Python: {get_python_version()}

INSTALLATION:

1. Install the wheel:

   pip install dist/akantu-{version}-*.whl

   Or with --user flag:

   pip install --user dist/akantu-{version}-*.whl

2. Verify installation:

   python3 -c "import akantu; print('Akantu successfully imported!')"

USAGE:

3. In your Python scripts:

   import akantu
   import numpy as np

   # Your Akantu code here...

NOTES:

- This wheel contains pre-built binaries for macOS
- Built with the dependencies that were present during compilation
- If you get import errors, you may need to install system dependencies:
  * Boost
  * Eigen
  * MUMPS (if parallel features used)
  * Scotch (if parallel features used)

TROUBLESHOOTING:

- ImportError: Make sure you're using Python {sys.version_info.major}.{sys.version_info.minor}
- Library errors: Check that system dependencies are available
- To uninstall: pip uninstall akantu

TRANSFERRING TO OTHER MACS:

This wheel should work on other Macs with:
- Same architecture ({platform.machine()})
- Same or compatible macOS version
- Same Python version ({sys.version_info.major}.{sys.version_info.minor})
- Similar system dependencies installed

Copy the wheel file to another Mac:
  scp dist/akantu-*.whl user@othermac:
  ssh user@othermac
  pip install --user akantu-*.whl

================================================================================
'''

    guide_path = DIST_DIR / "INSTALLATION_GUIDE.txt"
    guide_path.write_text(guide)
    print(f"\n✓ Created installation guide: {guide_path}")

def main():
    print("\n" + "=" * 60)
    print("Akantu Wheel Creator for macOS")
    print("=" * 60)

    # Verify paths exist
    if not AKANTU_ROOT.exists():
        print(f"✗ Error: Akantu root not found at {AKANTU_ROOT}")
        return 1

    if not PYTHON_BUILD_DIR.exists():
        print(f"✗ Error: Python build directory not found at {PYTHON_BUILD_DIR}")
        print("  Make sure Akantu is built with Python interface enabled")
        return 1

    print(f"✓ Akantu root: {AKANTU_ROOT}")
    print(f"✓ Build directory: {BUILD_DIR}")
    print(f"✓ Python module: {PYTHON_BUILD_DIR}")

    # Create wheel structure
    if not create_wheel_structure():
        return 1

    # Create setup.py
    version = create_setup_py()

    # Build wheel
    if not build_wheel():
        return 1

    # Create installation guide
    create_installation_guide(version)

    print("\n" + "=" * 60)
    print("SUCCESS!")
    print("=" * 60)
    print(f"\nWheel location: {DIST_DIR}")
    print(f"Installation guide: {DIST_DIR}/INSTALLATION_GUIDE.txt")
    print("\nTo install:")
    print(f"  pip install {DIST_DIR}/akantu-*.whl")
    print("\nTo verify:")
    print("  python3 -c 'import akantu; print(\"Success!\")'")
    print()

    return 0

if __name__ == "__main__":
    sys.exit(main())
