#!/bin/bash
# Setup environment for Akantu on macOS Apple Silicon
# Source this file before using akantu: source setup_akantu_env.sh

# Set library paths for Scotch error libraries (required for akantu to work)
export DYLD_INSERT_LIBRARIES="/opt/homebrew/lib/libscotcherr.dylib:/opt/homebrew/lib/libscotcherrexit.dylib"
export DYLD_LIBRARY_PATH="/opt/homebrew/lib:/usr/local/lib:${DYLD_LIBRARY_PATH}"

# Set Python path to find akantu module
export PYTHONPATH="/usr/local/lib/python3.10/site-packages:${PYTHONPATH}"

echo "✓ Akantu environment configured!"
echo "  - DYLD_INSERT_LIBRARIES set for Scotch error libraries"
echo "  - DYLD_LIBRARY_PATH includes /opt/homebrew/lib and /usr/local/lib"
echo "  - PYTHONPATH includes /usr/local/lib/python3.10/site-packages"
echo ""
echo "You can now use akantu with: python -c 'import akantu'"
echo "Make sure to use conda python (/Users/allamaprabhuani/miniconda3/bin/python)"
