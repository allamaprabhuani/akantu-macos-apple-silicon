#!/bin/bash
set -e

echo "=== Akantu GitLab Setup Script ==="
echo ""

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the current directory
REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

echo -e "${YELLOW}Step 1: Enter your GitLab repository URL${NC}"
echo "Example: https://gitlab.com/your-username/akantu-macos-apple-silicon.git"
read -p "GitLab URL: " GITLAB_URL

if [ -z "$GITLAB_URL" ]; then
    echo "Error: GitLab URL cannot be empty"
    exit 1
fi

echo ""
echo -e "${YELLOW}Step 2: Creating local branch for your work${NC}"

# Create a new branch from current detached HEAD
BRANCH_NAME="macos-apple-silicon-build"
git checkout -b "$BRANCH_NAME" 2>/dev/null || git checkout "$BRANCH_NAME"

echo -e "${GREEN}[OK] Created/checked out branch: $BRANCH_NAME${NC}"

echo ""
echo -e "${YELLOW}Step 3: Staging changes${NC}"

# Add documentation files
git add BUILD_MACOS_README.md CREATE_BACKUP.md .gitignore

# Add macOS build fixes
git add cmake/Modules/FindMumps.cmake
git add src/model/model.hh

echo -e "${GREEN}[OK] Staged macOS build fixes${NC}"

# Add phase field modifications
git add python/py_phasefield.cc
git add src/model/phase_field/

echo -e "${GREEN}[OK] Staged phase field modifications${NC}"

echo ""
echo -e "${YELLOW}Step 4: Creating commit${NC}"

# Create a comprehensive commit message
git commit -m "macOS Apple Silicon build support with phase field enhancements

This commit adds complete support for building Akantu on macOS with Apple Silicon,
along with tension-compression asymmetry features for phase field models.

macOS Build Fixes (Required):
- cmake/Modules/FindMumps.cmake: Fixed MUMPS detection for homebrew's MPI-enabled build
  * Added macOS-specific workaround for runtime test failures
  * Added automatic dependency resolution for BLAS, LAPACK, ScaLAPACK, gfortran
- src/model/model.hh: Fixed conditional compilation for CouplerSolidCohesiveContactOptions
  * Added AKANTU_COHESIVE_ELEMENT guard to prevent compilation errors

Phase Field Features (Optional):
- Added phi_minus and phi_plus fields for tension-compression asymmetry
- Implemented computePhiMinusOnQuad() in energy split classes
- Extended Python bindings for accessing new fields
- Modified quadratic phase field to compute asymmetric energy contributions

Documentation:
- BUILD_MACOS_README.md: Complete build guide with troubleshooting
- CREATE_BACKUP.md: Backup and restore strategies
- .gitignore: Added macOS-specific ignore patterns

Tested on:
- macOS 26.0.1 (Build 25A362)
- Architecture: arm64 (Apple Silicon)
- Compiler: GCC 15.2.0 (Homebrew)
- Python: 3.10.18
- Branch: origin/features/52-anisotropic-and-at1-phase-field-models
- Commit: 8c5e82d86" || echo "Commit already exists or no changes to commit"

echo -e "${GREEN}[OK] Committed changes${NC}"

echo ""
echo -e "${YELLOW}Step 5: Adding GitLab remote${NC}"

# Add remote (suppress error if already exists)
git remote add mygitlab "$GITLAB_URL" 2>/dev/null || {
    echo "Remote 'mygitlab' already exists, updating URL..."
    git remote set-url mygitlab "$GITLAB_URL"
}

echo -e "${GREEN}[OK] Added/updated remote: mygitlab${NC}"

echo ""
echo -e "${YELLOW}Step 6: Ready to push!${NC}"
echo ""
echo "Your local branch '$BRANCH_NAME' is ready to push to GitLab."
echo ""
echo "To push to GitLab, run:"
echo "  git push -u mygitlab $BRANCH_NAME"
echo ""
echo "After pushing, you can:"
echo "1. View your changes on GitLab"
echo "2. Create additional branches if needed"
echo "3. Set up CI/CD pipelines"
echo ""
echo -e "${GREEN}[OK] Setup complete!${NC}"
echo ""
echo "Current remotes:"
git remote -v
