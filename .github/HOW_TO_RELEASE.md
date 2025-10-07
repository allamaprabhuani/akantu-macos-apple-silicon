# How to Release Akantu Wheel on GitHub

## Prerequisites

- GitHub CLI: `brew install gh`
- Wheel file built: `packaging/create_standalone_wheel.py`
- Changes committed to repository

## Quick Release Process

```bash
# 1. Build wheel
python3 packaging/create_standalone_wheel.py

# 2. Get checksum
shasum -a 256 dist/akantu-5.0.7-*.whl

# 3. Create release
gh release create v5.0.7-macos-wheel \
  --title "Akantu v5.0.7 - Pre-built Wheel (Python 3.10)" \
  --notes-file docs/RELEASE_NOTES.md \
  --target macos-apple-silicon-build \
  dist/akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl
```

## Release Notes Template

Create `docs/RELEASE_NOTES.md`:

```markdown
# Akantu v5.0.7 - Pre-built Wheel for macOS Apple Silicon

## Quick Install

pip install --user akantu-5.0.7-cp310-none-macosx_26_0_arm64.whl

## Requirements

- macOS 11.0+, Apple Silicon
- Python 3.10
- Homebrew dependencies (see docs)

## SHA256

[checksum from step 2]

## Documentation

- Installation: docs/QUICK_START_WHEEL.md
- Technical: docs/README_WHEEL.md
```

## Update Release

```bash
# Delete old release
gh release delete v5.0.7-macos-wheel -y

# Create new one
gh release create v5.0.7-macos-wheel [...]
```

## For Different Python Versions

Build multiple wheels:

```bash
python3.10 packaging/create_standalone_wheel.py
python3.11 packaging/create_standalone_wheel.py
python3.12 packaging/create_standalone_wheel.py

# Upload all
gh release upload v5.0.7-macos-wheel dist/*.whl
```

## Checklist

- [ ] Wheel built and tested
- [ ] Checksum calculated and documented
- [ ] Release notes prepared
- [ ] Repository changes committed
- [ ] Release created on GitHub
- [ ] Download link works
- [ ] README updated with release link
