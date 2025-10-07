#!/usr/bin/env python3
"""
Create a standalone wheel from Akantu build with all dependencies bundled

This script:
1. Copies all required dynamic libraries
2. Fixes library paths using install_name_tool
3. Creates a completely standalone wheel

Usage:
    python3 create_standalone_wheel.py
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
WHEEL_BUILD_DIR = Path(__file__).parent / "wheel_build_standalone"
DIST_DIR = Path(__file__).parent / "dist"

def run_command(cmd, capture=True):
    """Run a shell command"""
    if capture:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.stdout.strip()
    else:
        subprocess.run(cmd, shell=True, check=True)

def get_akantu_version():
    """Extract version from CMakeLists.txt"""
    cmake_file = AKANTU_ROOT / "CMakeLists.txt"
    if cmake_file.exists():
        content = cmake_file.read_text()
        match = re.search(r'project\([^)]*VERSION\s+([0-9.]+)', content, re.IGNORECASE)
        if match:
            return match.group(1)
    return "5.0.7"

def get_dependencies(lib_path):
    """Get all dylib dependencies using otool"""
    output = run_command(f'otool -L "{lib_path}"')
    deps = []
    for line in output.split('\n')[1:]:  # Skip first line (the file itself)
        line = line.strip()
        if not line:
            continue
        # Extract library path (before the first parenthesis)
        lib = line.split('(')[0].strip()
        # Only include non-system libraries
        if lib.startswith('/opt/homebrew') or lib.startswith('@rpath') or lib.startswith('@loader_path'):
            deps.append(lib)
    return deps

def resolve_library_path(lib_name, search_paths):
    """Resolve @rpath or library name to actual path"""
    if lib_name.startswith('@rpath/'):
        base_name = lib_name.replace('@rpath/', '')
        for search_path in search_paths:
            full_path = Path(search_path) / base_name
            if full_path.exists():
                return str(full_path)
    elif lib_name.startswith('@loader_path/'):
        # Relative to the loading library - we'll handle this differently
        return None
    elif os.path.exists(lib_name):
        return lib_name
    return None

def collect_all_dependencies(initial_lib, search_paths):
    """Recursively collect all dependencies"""
    all_libs = {}
    to_process = [initial_lib]
    processed = set()

    while to_process:
        lib = to_process.pop(0)
        if lib in processed:
            continue

        processed.add(lib)

        # Resolve the actual path
        lib_path = resolve_library_path(lib, search_paths)
        if not lib_path or not os.path.exists(lib_path):
            continue

        # Get the library name without version
        lib_name = Path(lib_path).name
        all_libs[lib_name] = lib_path

        # Get this library's dependencies
        deps = get_dependencies(lib_path)
        for dep in deps:
            if dep not in processed:
                to_process.append(dep)

    return all_libs

def fix_library_paths(lib_path, libs_dir):
    """Fix library paths to use @loader_path"""
    print(f"  Fixing paths in {Path(lib_path).name}...")

    # Get current dependencies
    deps = get_dependencies(lib_path)

    for dep in deps:
        dep_basename = Path(dep.replace('@rpath/', '')).name

        # Change the reference to use @loader_path
        new_path = f"@loader_path/.dylibs/{dep_basename}"

        try:
            run_command(
                f'install_name_tool -change "{dep}" "{new_path}" "{lib_path}"',
                capture=False
            )
        except:
            # Try with the full path if @rpath didn't work
            if '/opt/homebrew' in dep:
                try:
                    run_command(
                        f'install_name_tool -change "{dep}" "{new_path}" "{lib_path}"',
                        capture=False
                    )
                except:
                    pass

    # Fix the library's own ID if it's a dylib
    if str(lib_path).endswith('.dylib'):
        lib_name = Path(lib_path).name
        run_command(
            f'install_name_tool -id "@loader_path/.dylibs/{lib_name}" "{lib_path}"',
            capture=False
        )

def create_standalone_wheel():
    """Create standalone wheel with all dependencies bundled"""
    print("=" * 70)
    print("Creating Standalone Akantu Wheel for macOS")
    print("=" * 70)

    # Clean and create directories
    if WHEEL_BUILD_DIR.exists():
        shutil.rmtree(WHEEL_BUILD_DIR)
    WHEEL_BUILD_DIR.mkdir(parents=True)

    if DIST_DIR.exists():
        for f in DIST_DIR.glob("akantu-*.whl"):
            f.unlink()
    else:
        DIST_DIR.mkdir(parents=True)

    # Create package directory
    package_dir = WHEEL_BUILD_DIR / "akantu"
    package_dir.mkdir()

    print("\n1. Copying Python module files...")
    # Copy __init__.py and modify it to preload scotch libraries
    init_content = (PYTHON_BUILD_DIR / "__init__.py").read_text()

    # Add preload code at the beginning
    preload_code = """# Preload Scotch error libraries for symbol resolution
import ctypes
import os
import sys
_lib_dir = os.path.join(os.path.dirname(__file__), '.dylibs')
if os.path.exists(_lib_dir):
    # Pre-load scotch error libraries with RTLD_GLOBAL to resolve symbols
    for _lib in ['libscotcherr.7.0.dylib', 'libscotcherrexit.7.0.dylib']:
        _lib_path = os.path.join(_lib_dir, _lib)
        if os.path.exists(_lib_path):
            try:
                ctypes.CDLL(_lib_path, mode=ctypes.RTLD_GLOBAL)
            except:
                pass  # Silently ignore if preload fails

"""

    # Insert after the copyright/license header
    lines = init_content.split('\n')
    insert_pos = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.strip().startswith('#') and not line.strip().startswith('"""') and not line.strip().startswith("'''"):
            insert_pos = i
            break

    lines.insert(insert_pos, preload_code)
    modified_init = '\n'.join(lines)

    (package_dir / "__init__.py").write_text(modified_init)
    print(f"  ✓ Created modified __init__.py with scotch preload")

    # Copy the compiled extension
    so_files = list(PYTHON_BUILD_DIR.glob("*.so"))
    main_so = None
    for so_file in so_files:
        main_so = package_dir / so_file.name
        shutil.copy2(so_file, main_so)
        print(f"  ✓ Copied {so_file.name}")

    print("\n2. Collecting all dynamic library dependencies...")
    # Search paths for @rpath resolution
    search_paths = [
        BUILD_DIR / "src",
        BUILD_DIR / "third-party" / "iohelper" / "src",
        "/opt/homebrew/lib",
        "/opt/homebrew/opt/scotch/lib",
        "/opt/homebrew/opt/open-mpi/lib",
        "/opt/homebrew/opt/brewsci-mumps/lib",
        "/opt/homebrew/opt/gcc/lib/gcc/current",
        "/opt/homebrew/opt/scalapack/lib",
        "/opt/homebrew/opt/openblas/lib",
    ]
    search_paths = [str(p) for p in search_paths if Path(p).exists()]

    # Collect all dependencies recursively
    all_libs = collect_all_dependencies(str(main_so), search_paths)

    # Add all scotch libraries explicitly (they're loaded at runtime)
    scotch_libs = [
        "/opt/homebrew/lib/libscotch.7.0.dylib",
        "/opt/homebrew/lib/libscotcherr.7.0.dylib",
        "/opt/homebrew/lib/libscotcherrexit.7.0.dylib",
    ]
    for scotch_lib in scotch_libs:
        if Path(scotch_lib).exists():
            lib_name = Path(scotch_lib).name
            if lib_name not in all_libs:
                all_libs[lib_name] = scotch_lib
                print(f"  + Added runtime dependency: {lib_name}")

    print(f"  Found {len(all_libs)} unique libraries")

    # Create .dylibs directory
    libs_dir = package_dir / ".dylibs"
    libs_dir.mkdir()

    print("\n3. Copying libraries...")
    for lib_name, lib_path in all_libs.items():
        dest = libs_dir / lib_name
        shutil.copy2(lib_path, dest)
        # Make it writable so we can modify it
        os.chmod(dest, 0o755)
        print(f"  ✓ {lib_name}")

    print("\n4. Fixing library paths...")
    # Fix paths in the main .so file
    fix_library_paths(main_so, libs_dir)

    # Fix paths in all copied libraries
    for lib_file in libs_dir.glob("*.dylib"):
        fix_library_paths(lib_file, libs_dir)

    print("\n5. Creating setup.py...")
    version = get_akantu_version()

    import platform
    mac_ver = platform.mac_ver()[0].split('.')
    arch = platform.machine()
    if arch == 'arm64':
        plat_name = f"macosx_{mac_ver[0]}_0_arm64"
    else:
        plat_name = f"macosx_{mac_ver[0]}_0_x86_64"

    setup_content = f'''#!/usr/bin/env python
from setuptools import setup, find_packages
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel

class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        _bdist_wheel.finalize_options(self)
        self.root_is_pure = False

    def get_tag(self):
        python = "py3"
        abi = "none"
        plat = "{plat_name}"
        return python, abi, plat

setup(
    name="akantu",
    version="{version}",
    url="https://akantu.ch",
    author="Nicolas Richart",
    author_email="nicolas.richart@epfl.ch",
    description="Akantu finite element library for macOS Apple Silicon",
    long_description="Pre-compiled Akantu library for macOS with dependencies bundled.",
    platforms="macOS",
    license="L-GPLv3",
    install_requires=["numpy", "scipy"],
    packages=find_packages(),
    package_data={{
        "akantu": ["*.so", ".dylibs/*"]
    }},
    include_package_data=True,
    zip_safe=False,
    cmdclass={{"bdist_wheel": bdist_wheel}},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
    ],
)
'''
    (WHEEL_BUILD_DIR / "setup.py").write_text(setup_content)
    print(f"  ✓ Created setup.py (version {version}, platform {plat_name})")

    print("\n6. Building wheel...")
    os.chdir(WHEEL_BUILD_DIR)
    result = subprocess.run(
        [sys.executable, "setup.py", "bdist_wheel"],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("Error building wheel:")
        print(result.stderr)
        return False

    # Move wheel to dist
    wheel_dist = WHEEL_BUILD_DIR / "dist"
    for wheel_file in wheel_dist.glob("*.whl"):
        dest = DIST_DIR / wheel_file.name
        if dest.exists():
            dest.unlink()
        shutil.copy2(wheel_file, dest)
        print(f"\n✓ Wheel created: {dest}")
        print(f"  Size: {dest.stat().st_size / (1024*1024):.1f} MB")
        print(f"  Libraries bundled: {len(all_libs)}")

    return True

def main():
    print("\nThis script creates a standalone wheel with all dependencies bundled.\n")

    if not PYTHON_BUILD_DIR.exists():
        print(f"✗ Error: Python build directory not found at {PYTHON_BUILD_DIR}")
        return 1

    if not create_standalone_wheel():
        return 1

    print("\n" + "=" * 70)
    print("SUCCESS!")
    print("=" * 70)
    print(f"\nStandalone wheel created in: {DIST_DIR}")
    print("\nThis wheel includes all dependencies and should work on any")
    print("Apple Silicon Mac with Python 3.10, without requiring Homebrew.")
    print("\nTo install:")
    print(f"  pip install {DIST_DIR}/akantu-*.whl")
    print("\nTo test:")
    print("  python3 -c 'import akantu; print(\"Success!\")'")
    print()

    return 0

if __name__ == "__main__":
    sys.exit(main())
