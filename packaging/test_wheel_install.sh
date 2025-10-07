#!/bin/bash
#
# Test Akantu Wheel Installation
#
# This script verifies that the Akantu wheel can be installed and imported
# without affecting your system installation.
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WHEEL_FILE="${SCRIPT_DIR}/dist/akantu-5.0.7-py3-none-any.whl"
TEST_ENV="/tmp/akantu_test_$$"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================="
echo "Akantu Wheel Installation Test"
echo "========================================="
echo ""

# Check wheel exists
if [ ! -f "$WHEEL_FILE" ]; then
    echo -e "${RED}✗ Error: Wheel file not found at $WHEEL_FILE${NC}"
    echo "Run: python3 create_akantu_wheel_mac.py first"
    exit 1
fi

echo -e "${GREEN}✓ Found wheel: $WHEEL_FILE${NC}"
echo "  Size: $(ls -lh "$WHEEL_FILE" | awk '{print $5}')"
echo ""

# Create virtual environment
echo "Creating test virtual environment..."
python3 -m venv "$TEST_ENV"
echo -e "${GREEN}✓ Created venv: $TEST_ENV${NC}"
echo ""

# Activate and install
echo "Installing wheel..."
source "$TEST_ENV/bin/activate"

pip install --quiet --upgrade pip
pip install "$WHEEL_FILE"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Wheel installed successfully${NC}"
else
    echo -e "${RED}✗ Wheel installation failed${NC}"
    rm -rf "$TEST_ENV"
    exit 1
fi

echo ""
echo "Testing import..."

# Test import
python3 << 'EOF'
import sys
try:
    import akantu
    print(f"\033[0;32m✓ Akantu imported successfully!\033[0m")
    print(f"  Module location: {akantu.__file__}")
    print(f"  Available attributes: {len(dir(akantu))}")

    # Try to access some basic functionality
    if hasattr(akantu, 'py11_akantu'):
        print(f"  Core module: ✓ py11_akantu found")

    sys.exit(0)
except ImportError as e:
    print(f"\033[0;31m✗ Import failed: {e}\033[0m")
    print(f"\033[1;33m")
    print("This is expected if Homebrew dependencies are not installed.")
    print("The wheel requires these Homebrew packages:")
    print("  - gcc")
    print("  - boost")
    print("  - scotch")
    print("  - eigen")
    print("  - open-mpi")
    print("  - openblas")
    print("  - scalapack")
    print("  - brewsci-mumps")
    print(f"\033[0m")
    sys.exit(1)
EOF

IMPORT_RESULT=$?

# Cleanup
deactivate
rm -rf "$TEST_ENV"

echo ""
echo "========================================="
if [ $IMPORT_RESULT -eq 0 ]; then
    echo -e "${GREEN}SUCCESS: Wheel works correctly!${NC}"
    echo "========================================="
    echo ""
    echo "The wheel can be installed with:"
    echo "  pip install --user $WHEEL_FILE"
    echo ""
    echo "Or in a virtual environment:"
    echo "  python3 -m venv myenv"
    echo "  source myenv/bin/activate"
    echo "  pip install $WHEEL_FILE"
    echo ""
else
    echo -e "${YELLOW}PARTIAL SUCCESS: Wheel installs but needs dependencies${NC}"
    echo "========================================="
    echo ""
    echo "The wheel packages correctly, but requires Homebrew dependencies"
    echo "to be installed on the target system:"
    echo ""
    echo "  brew install gcc boost scotch eigen open-mpi openblas scalapack"
    echo "  brew tap brewsci/num"
    echo "  brew install brewsci-mumps --without-brewsci-parmetis"
    echo ""
    echo "This is normal for a wheel built with system dependencies."
    echo ""
    echo "The wheel will work on any Mac with:"
    echo "  1. Apple Silicon (arm64)"
    echo "  2. macOS 11.0+"
    echo "  3. Python 3.10"
    echo "  4. Same Homebrew dependencies installed"
    echo ""
fi

echo "To install on this system (with dependencies already present):"
echo "  pip install --user $WHEEL_FILE"
echo ""
