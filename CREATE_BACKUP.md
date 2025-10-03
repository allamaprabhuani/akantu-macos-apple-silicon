# Creating Portable Backups of Akantu Installation

This guide explains how to create backups of your working Akantu installation for easy restoration on the same or similar systems.

## Option 1: Create a Tarball of Installed Files (Simplest)

Package all installed files into a single archive:

```bash
# Create backup directory
mkdir -p ~/akantu_backups

# Create tarball of installed files
cd /usr/local
sudo tar -czf ~/akantu_backups/akantu_install_$(date +%Y%m%d).tar.gz \
    lib/libakantu* \
    lib/libiohelper* \
    lib/python3.10/site-packages/akantu \
    include/akantu* \
    bin/*akantu* 2>/dev/null || true

# Also backup the source with modifications
cd /Users/allamaprabhuani/akantu_mac_oct
tar -czf ~/akantu_backups/akantu_source_$(date +%Y%m%d).tar.gz \
    --exclude='akantu/build' \
    --exclude='akantu/.git' \
    akantu/

echo "Backups created in ~/akantu_backups/"
ls -lh ~/akantu_backups/
```

### Restore from Tarball

```bash
# Restore installed files
cd /usr/local
sudo tar -xzf ~/akantu_backups/akantu_install_YYYYMMDD.tar.gz

# Fix library install names (if needed)
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib

# Add to PYTHONPATH in ~/.zshrc
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"
```

## Option 2: Create a Relocatable Build (Advanced)

Build Akantu with relocatable binaries that can be moved anywhere:

```bash
cd /Users/allamaprabhuani/akantu_mac_oct/akantu/build

# Reconfigure with relative RPATH
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/akantu_portable \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -DCMAKE_INSTALL_RPATH="@loader_path/../lib;@loader_path" \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=ON \
  -DMUMPS_INCLUDE_DIR=/opt/homebrew/opt/brewsci-mumps/include \
  -DMUMPS_LIBRARY_DMUMPS=/opt/homebrew/opt/brewsci-mumps/lib/libdmumps.dylib \
  -DMUMPS_LIBRARY_COMMON=/opt/homebrew/opt/brewsci-mumps/lib/libmumps_common.dylib \
  -DMUMPS_LIBRARY_PORD=/opt/homebrew/opt/brewsci-mumps/lib/libpord.dylib \
  -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"

make -j$(sysctl -n hw.ncpu)
make install

# Create portable package
cd $HOME
tar -czf ~/akantu_backups/akantu_portable_$(date +%Y%m%d).tar.gz akantu_portable/
```

### Use Portable Installation

```bash
# Extract anywhere
tar -xzf akantu_portable_YYYYMMDD.tar.gz -C ~/wherever/

# Add to PYTHONPATH
export PYTHONPATH="$HOME/wherever/akantu_portable/lib/python3.10/site-packages:$PYTHONPATH"
```

**Note**: This still requires homebrew dependencies (MUMPS, Scotch, etc.) to be installed.

## Option 3: Create a Conda Package (Best for Distribution)

Create a conda package that includes all dependencies:

```bash
# Install conda-build
conda install conda-build

# Create conda recipe directory
mkdir -p ~/akantu_conda_recipe
cd ~/akantu_conda_recipe

# Create meta.yaml (recipe file)
cat > meta.yaml << 'EOF'
package:
  name: akantu
  version: "5.0.7"

build:
  number: 0
  skip: True  # [not osx]

requirements:
  build:
    - {{ compiler('c') }}
    - {{ compiler('cxx') }}
    - {{ compiler('fortran') }}
    - cmake
    - make
  host:
    - python
    - numpy
    - scipy
    - boost
    - eigen
    - open-mpi
    - openblas
    - scalapack
  run:
    - python
    - numpy
    - scipy
    - open-mpi
    - openblas

about:
  home: https://gitlab.com/akantu/akantu
  license: LGPL-3.0
  summary: Swiss-Made Open-Source Finite-Element Library

EOF

# Create build script
cat > build.sh << 'EOF'
#!/bin/bash

mkdir build && cd build

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=OFF \
  -DAKANTU_PYTHON_INTERFACE=ON

make -j${CPU_COUNT}
make install
EOF

chmod +x build.sh

# Build conda package
conda build . --output-folder ~/akantu_backups/conda_packages/
```

### Install Conda Package

```bash
# Install from local package
conda install --use-local akantu

# Or create a local channel
conda index ~/akantu_backups/conda_packages/
conda install -c file:///$HOME/akantu_backups/conda_packages akantu
```

## Option 4: Docker Image (Most Portable)

Create a Docker image with full working environment:

```bash
# Create Dockerfile
cat > ~/akantu_docker/Dockerfile << 'EOF'
FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    cmake build-essential \
    libboost-dev libeigen3-dev \
    mpi-default-dev libmumps-dev libscotch-dev \
    python3 python3-pip python3-numpy python3-scipy

# Copy source code
COPY akantu /akantu

# Build Akantu
WORKDIR /akantu/build
RUN cmake .. && make -j$(nproc) && make install

# Set up Python environment
ENV PYTHONPATH=/usr/local/lib/python3.10/site-packages:$PYTHONPATH

WORKDIR /workspace
CMD ["/bin/bash"]
EOF

# Build Docker image
docker build -t akantu:macos-apple-silicon ~/akantu_docker/

# Save Docker image
docker save akantu:macos-apple-silicon | gzip > ~/akantu_backups/akantu_docker_$(date +%Y%m%d).tar.gz

# Load on another system
docker load < akantu_docker_YYYYMMDD.tar.gz
```

## Option 5: Git Repository + Build Script (For Development)

Keep source + automated build script:

```bash
# Create installation script
cat > ~/akantu_mac_oct/akantu/install_akantu.sh << 'EOF'
#!/bin/bash
set -e

echo "Installing Akantu on macOS..."

# Check dependencies
echo "Checking homebrew dependencies..."
brew list gcc boost scotch eigen open-mpi openblas scalapack brewsci-mumps || {
    echo "Installing missing dependencies..."
    brew install gcc boost scotch eigen open-mpi openblas scalapack
    brew tap brewsci/num
    brew install brewsci-mumps --without-brewsci-parmetis
}

# Build
cd "$(dirname "$0")"
mkdir -p build && cd build

echo "Configuring CMake..."
cmake .. \
  -DAKANTU_PARALLEL=ON \
  -DAKANTU_USE_SYSTEM_MUMPS=ON \
  -DMUMPS_INCLUDE_DIR=/opt/homebrew/opt/brewsci-mumps/include \
  -DMUMPS_LIBRARY_DMUMPS=/opt/homebrew/opt/brewsci-mumps/lib/libdmumps.dylib \
  -DMUMPS_LIBRARY_COMMON=/opt/homebrew/opt/brewsci-mumps/lib/libmumps_common.dylib \
  -DMUMPS_LIBRARY_PORD=/opt/homebrew/opt/brewsci-mumps/lib/libpord.dylib \
  -DSCOTCH_LIBRARY="/opt/homebrew/lib/libscotch.dylib;/opt/homebrew/lib/libscotcherr.dylib;/opt/homebrew/lib/libscotcherrexit.dylib"

echo "Building..."
make -j$(sysctl -n hw.ncpu)

echo "Installing (requires sudo)..."
sudo make install

echo "Fixing library install names..."
sudo install_name_tool -id /usr/local/lib/libakantu.5.0.dylib /usr/local/lib/libakantu.5.0.dylib
sudo install_name_tool -id /usr/local/lib/libiohelper.dylib /usr/local/lib/libiohelper.dylib

echo "Fixing Python module..."
sudo cp /usr/local/lib/python3.10/site-packages/akantu/akantu/__init__.py \
        /usr/local/lib/python3.10/site-packages/akantu/__init__.py 2>/dev/null || true

# Fix mpi4py import
sudo sed -i '' 's/except Exception:/except (Exception, ModuleNotFoundError):/' \
    /usr/local/lib/python3.10/site-packages/akantu/__init__.py 2>/dev/null || true

echo "Setting up environment..."
if ! grep -q "PYTHONPATH.*python3.10/site-packages" ~/.zshrc; then
    echo 'export PYTHONPATH="/usr/local/lib/python3.10/site-packages:$PYTHONPATH"' >> ~/.zshrc
    echo "Added PYTHONPATH to ~/.zshrc"
fi

echo ""
echo "Installation complete!"
echo "Please run: source ~/.zshrc"
echo "Then test with: python -c 'import akantu; print(len(dir(akantu)))'"
EOF

chmod +x ~/akantu_mac_oct/akantu/install_akantu.sh
```

### Create Full Backup

```bash
# Commit all changes
cd ~/akantu_mac_oct/akantu
git add -A
git commit -m "macOS build with phase field modifications"

# Create bundle (includes all git history)
git bundle create ~/akantu_backups/akantu_macos_$(date +%Y%m%d).bundle --all

# Or push to a private repository
git remote add backup <your-private-repo-url>
git push backup --all
```

### Restore from Git Bundle

```bash
# Clone from bundle
git clone ~/akantu_backups/akantu_macos_YYYYMMDD.bundle akantu

# Run installation script
cd akantu
./install_akantu.sh
```

## Recommended Approach

**For your use case**, I recommend **Option 5 (Git + Install Script)** combined with **Option 1 (Tarball backup)**:

1. Keep source in git with modifications
2. Use the install script for easy rebuilds
3. Keep a tarball of working installation for quick restore

```bash
# Quick backup command
mkdir -p ~/akantu_backups
cd /Users/allamaprabhuani/akantu_mac_oct/akantu
git bundle create ~/akantu_backups/akantu_source_$(date +%Y%m%d).bundle --all
cd /usr/local && sudo tar -czf ~/akantu_backups/akantu_binaries_$(date +%Y%m%d).tar.gz \
    lib/libakantu* lib/libiohelper* lib/python3.10/site-packages/akantu
```

This gives you:
- Source code with all modifications (git bundle)
- Compiled binaries for instant restore (tarball)
- Easy rebuild capability (install script)

## Restore Checklist

Whenever you restore on a fresh system:

1. [x] Install homebrew dependencies (gcc, boost, mumps, etc.)
2. [x] Restore source or extract tarball
3. [x] Run install script or extract binaries
4. [x] Fix library install names
5. [x] Add PYTHONPATH to shell config
6. [x] Remove mpi4py if present
7. [x] Test: `python -c "import akantu; print(len(dir(akantu)))"`
