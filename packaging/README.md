# Akantu Wheel Packaging

This directory contains scripts for creating Python wheel distributions of Akantu for macOS Apple Silicon.

## Quick Usage

```bash
python3 create_standalone_wheel.py
```

The wheel will be created in `../dist/`

## Scripts

- **create_standalone_wheel.py** - Creates wheel with all dependencies bundled (recommended)
- **create_akantu_wheel_mac.py** - Creates minimal wheel requiring system dependencies
- **test_wheel_install.sh** - Tests wheel installation in clean environment
- **CHECKSUM.txt** - SHA256 verification for distributed wheel

## Documentation

See [../docs/README_WHEEL.md](../docs/README_WHEEL.md) for complete information.
