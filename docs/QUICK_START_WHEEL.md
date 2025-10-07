# Quick Start: Akantu Wheel Installation

## What You Have

A Python wheel package of your Akantu build that can be installed in seconds instead of rebuilding for 1-2 hours.

**Wheel file**: `dist/akantu-5.0.7-py3-none-any.whl` (38.6 MB)

## Install on This Mac

```bash
pip install --user dist/akantu-5.0.7-py3-none-any.whl
```

## Verify It Works

```bash
python3 -c "import akantu; print('Akantu ready!')"
```

##  Usage

```python
import akantu

# Your code here
mesh = akantu.Mesh(spatial_dimension=2)
# ... rest of your simulation
```

## Use in a Virtual Environment

```bash
# Create environment
python3 -m venv my_akantu_env
source my_akantu_env/bin/activate

# Install
pip install /path/to/dist/akantu-5.0.7-py3-none-any.whl

# Use
python your_simulation.py
```

## Transfer to Another Mac

1. **Copy the wheel**:
   ```bash
   scp dist/akantu-5.0.7-py3-none-any.whl user@other-mac:
   ```

2. **On the other Mac**, install Homebrew dependencies:
   ```bash
   brew install gcc boost scotch eigen open-mpi openblas scalapack
   brew tap brewsci/num
   brew install brewsci-mumps --without-brewsci-parmetis
   ```

3. **Install the wheel**:
   ```bash
   pip install --user akantu-5.0.7-py3-none-any.whl
   ```

## Rebuild the Wheel After Updates

```bash
# After rebuilding Akantu
python3 create_standalone_wheel.py
```

## What's Different from Source Build?

| | Source Build | Wheel Install |
|---|---|---|
| **Time** | 1-2 hours | 10 seconds |
| **Command** | `make install` | `pip install` |
| **Location** | System-wide | User or venv |
| **Uninstall** | Manual | `pip uninstall akantu` |
| **Update** | Rebuild everything | Reinstall wheel |

## Requirements

- macOS 11.0+ (Apple Silicon)
- Python 3.10
- Homebrew with these packages:
  - gcc, boost, scotch, eigen
  - open-mpi, openblas, scalapack
  - brewsci-mumps

## That's It!

You now have a fast, portable way to install your Akantu build. No more waiting hours for compilation!

For more details, see [README_WHEEL.md](README_WHEEL.md)
