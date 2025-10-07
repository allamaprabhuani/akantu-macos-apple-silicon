# Quick Start: Akantu Wheel Installation

## Install Pre-built Wheel

```bash
# Download
wget https://github.com/allamaprabhuani/akantu-macos-apple-silicon/releases/download/v5.0.7-macos-wheel/akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl

# Install dependencies
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis

# Install wheel
pip install --user akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl

# Verify
python3 -c "import akantu; print('Installation successful')"
```

## Requirements

- macOS 11.0+, Apple Silicon
- Python 3.10 (check: `python3 --version`)
- Homebrew dependencies (above)

## Virtual Environment

```bash
python3 -m venv akantu_env
source akantu_env/bin/activate
pip install akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl
```

## Build Custom Wheel

For different Python versions:

```bash
python3 packaging/create_standalone_wheel.py
# Output: dist/akantu-5.0.7-cp3XX-none-macosx_26_0_arm64.whl
```

## Troubleshooting

**ImportError: symbol not found**
- Install Homebrew dependencies (brew install commands above)

**Wrong Python version**
- Wheel requires Python 3.10
- Build custom wheel for your version

See [README_WHEEL.md](README_WHEEL.md) for technical details.
