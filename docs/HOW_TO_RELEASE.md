# How to Release the Akantu Wheel on GitHub

## Overview

This guide shows you how to publish the pre-built Akantu wheel to your GitHub repository so other users can download and use it without building from source.

## Prerequisites

- GitHub account with access to https://github.com/allamaprabhuani/akantu-macos-apple-silicon
- `gh` CLI tool (GitHub CLI) installed: `brew install gh`
- Or use GitHub web interface

## Files to Include in Release

1. **The Wheel File**:
   - `dist/akantu-5.0.7-py3-none-any.whl` (38.6 MB)

2. **Documentation** (add to repository):
   - `QUICK_START_WHEEL.md`
   - `README_WHEEL.md`
   - `create_standalone_wheel.py`
   - `create_akantu_wheel_mac.py`
   - `test_wheel_install.sh`

3. **Release Notes**:
   - Use `GITHUB_RELEASE.md` as template

## Step-by-Step Release Process

### Option 1: Using GitHub Web Interface (Easiest)

1. **Prepare the files in your repository**:
   ```bash
   cd /Users/allamaprabhuani/akantu_mac_oct/akantu

   # Copy wheel packaging scripts to repo
   cp ../create_standalone_wheel.py .
   cp ../create_akantu_wheel_mac.py .
   cp ../test_wheel_install.sh .
   cp ../QUICK_START_WHEEL.md .
   cp ../README_WHEEL.md .
   ```

2. **Commit and push to GitHub**:
   ```bash
   git add create_standalone_wheel.py create_akantu_wheel_mac.py test_wheel_install.sh
   git add QUICK_START_WHEEL.md README_WHEEL.md
   git commit -m "Add pre-built wheel packaging scripts and documentation"
   git push github main
   ```

3. **Create GitHub Release**:
   - Go to https://github.com/allamaprabhuani/akantu-macos-apple-silicon/releases
   - Click "Draft a new release"
   - Set tag: `v5.0.7-wheel`
   - Release title: `Akantu v5.0.7 - Pre-built Wheel for macOS Apple Silicon`
   - Description: Copy content from `GITHUB_RELEASE.md`
   - Upload `dist/akantu-5.0.7-py3-none-any.whl`
   - Click "Publish release"

### Option 2: Using GitHub CLI (Faster)

```bash
cd /Users/allamaprabhuani/akantu_mac_oct/akantu

# Login to GitHub (if not already)
gh auth login

# Create release with the wheel attached
gh release create v5.0.7-wheel \
  --title "Akantu v5.0.7 - Pre-built Wheel for macOS Apple Silicon" \
  --notes-file ../GITHUB_RELEASE.md \
  ../dist/akantu-5.0.7-py3-none-any.whl
```

## Update README.md

Add a section to your main [README.md](akantu/README.md) about the pre-built wheel:

```markdown
## 🚀 Quick Install - Pre-built Wheel (New!)

**Don't want to build from source?** Download the pre-built wheel:

```bash
# Download from releases
wget https://github.com/allamaprabhuani/akantu-macos-apple-silicon/releases/download/v5.0.7-wheel/akantu-5.0.7-py3-none-any.whl

# Install (requires Homebrew dependencies)
brew install gcc boost scotch eigen open-mpi openblas scalapack
brew tap brewsci/num
brew install brewsci-mumps --without-brewsci-parmetis

pip install --user akantu-5.0.7-py3-none-any.whl
```

**Installation time**: 10 seconds vs 1-2 hours building from source!

See [QUICK_START_WHEEL.md](QUICK_START_WHEEL.md) for details.
```

## Recommended Repository Structure

```
akantu-macos-apple-silicon/
├── README.md                      # Main readme (update with wheel info)
├── QUICK_START_WHEEL.md          # Quick wheel installation guide
├── README_WHEEL.md               # Detailed wheel documentation
├── create_standalone_wheel.py    # Script to create wheel
├── create_akantu_wheel_mac.py    # Alternative wheel script
├── test_wheel_install.sh         # Test script
├── src/                          # Akantu source code
├── python/                       # Python bindings
└── ... (rest of Akantu files)
```

## After Release

### Announce to Users

1. **Update main README.md** with prominent wheel download section
2. **Pin the release** on GitHub (make it easy to find)
3. **Create a discussion** in GitHub Discussions about the wheel
4. **Update any documentation** that mentions "building from source" to include wheel option

### Monitor Feedback

1. Watch for issues related to wheel installation
2. Check if users have problems with specific macOS versions
3. Consider creating wheels for other Python versions (3.11, 3.12) if requested

## Updating the Wheel

When you update Akantu:

1. **Rebuild Akantu**:
   ```bash
   cd akantu/build
   make -j$(sysctl -n hw.ncpu)
   ```

2. **Recreate wheel**:
   ```bash
   cd ..
   python3 create_standalone_wheel.py
   ```

3. **Create new release**:
   ```bash
   gh release create v5.0.8-wheel \
     --title "Akantu v5.0.8 - Updated Wheel" \
     --notes "Updated Akantu build with latest changes" \
     dist/akantu-5.0.8-py3-none-any.whl
   ```

## Important Notes

### File Size Considerations

- The wheel is 38.6 MB (compressed)
- GitHub has a 2 GB file size limit per release, so no issues
- Users with slow internet may take a minute to download

### Platform Specificity

This wheel only works for:
- macOS 11.0+
- Apple Silicon (arm64)
- Python 3.10

Consider mentioning this prominently in documentation.

### License Compliance

The wheel bundles several libraries:
- Make sure your LICENSE file mentions bundled components
- Include attribution for Scotch, MUMPS, GCC runtime, etc.
- LGPLv3 allows distribution but requires source code availability

### Security

- The wheel contains pre-compiled binaries
- Users may want to verify checksums
- Consider providing SHA256 hash in release notes:

```bash
shasum -a 256 dist/akantu-5.0.7-py3-none-any.whl
```

## Testing Before Release

Run these checks:

```bash
# 1. Test wheel installation
./test_wheel_install.sh

# 2. Check wheel metadata
unzip -l dist/akantu-5.0.7-py3-none-any.whl

# 3. Verify all libraries are bundled
unzip -q dist/akantu-5.0.7-py3-none-any.whl -d /tmp/test_wheel
ls -lh /tmp/test_wheel/akantu/.dylibs/

# 4. Test in fresh venv
python3 -m venv /tmp/test_akantu
source /tmp/test_akantu/bin/activate
pip install dist/akantu-5.0.7-py3-none-any.whl
python3 -c "import akantu"
```

## Example Release Checklist

- [ ] Wheel file created and tested
- [ ] Documentation files added to repository
- [ ] README.md updated with wheel information
- [ ] Release notes prepared (GITHUB_RELEASE.md)
- [ ] SHA256 checksum calculated
- [ ] Wheel uploaded to GitHub release
- [ ] Release published and pinned
- [ ] Users notified (discussions, social media, etc.)

## Questions?

If users have issues:
1. Direct them to QUICK_START_WHEEL.md
2. Check if Homebrew dependencies are installed
3. Verify Python version compatibility
4. Open issue on GitHub for tracking

---

**Ready to release?** Follow the steps above and your wheel will be available to the community!
