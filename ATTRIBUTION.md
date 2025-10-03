# Attribution and Modifications

## Acknowledgments

This work builds upon the **Akantu** finite element library, developed by the Computational Solid Mechanics Laboratory (LSMS) at École Polytechnique Fédérale de Lausanne (EPFL). I gratefully acknowledge the original developers for making this powerful tool available under an open-source license.

## Original Source

**Repository**: https://gitlab.com/akantu/akantu
**Branch**: `features/52-anisotropic-and-at1-phase-field-models`
**Commit**: 8c5e82d86
**License**: GNU Lesser General Public License v3.0 (LGPLv3)
**Copyright**: (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
**Laboratory**: LSMS - Laboratoire de Simulation en Mécanique des Solides

### Citation

When using Akantu in academic work, please cite:

```bibtex
@article{akantu2024,
  title={Akantu: an HPC finite-element library for contact and dynamic fracture simulations},
  author={Richart, Nicolas and Molinari, Jean-Fran{\c{c}}ois and others},
  journal={The Journal of Open Source Software},
  year={2024},
  publisher={The Open Journal}
}
```

Documentation: https://akantu.readthedocs.io

## Modifications and Contributions

This fork contains two categories of modifications: (1) platform-specific build fixes required for macOS Apple Silicon compatibility, and (2) research-oriented enhancements to the phase field fracture framework.

### 1. Platform-Specific Build Adaptations (Technical Requirement)

The following adaptations were necessary to enable compilation on macOS Apple Silicon (arm64) architecture. These modifications address platform-specific differences in library packaging and compiler behavior between Linux and macOS environments.

#### `cmake/Modules/FindMumps.cmake` (+17 lines)

**Technical Issue**: The Homebrew distribution of MUMPS on macOS is compiled with MPI support rather than as a sequential library. During CMake configuration, the test binary compiles successfully but fails at runtime due to MPI initialization requirements in the testing environment, causing the configuration process to abort after multiple retries.

**Modifications Implemented**:
1. Lines ~282-288: Added macOS-specific logic to accept successful compilation as sufficient evidence of MUMPS availability, bypassing the runtime test requirement after 5 failed attempts
2. Lines ~361-372: Extended automatic dependency resolution to explicitly link required Fortran and numerical libraries (mumps_common, pord, BLAS, LAPACK, ScaLAPACK, gfortran) that are not automatically detected in the Homebrew package structure

**Rationale**: This approach maintains compatibility with the original Linux build system while accommodating macOS-specific library packaging conventions.

#### `src/model/model.hh` (+2 lines)

**Technical Issue**: Conditional compilation error arising from inconsistent preprocessor guard usage. The `CouplerSolidCohesiveContactOptions` class is only defined when the `AKANTU_COHESIVE_ELEMENT` feature is enabled, but the code attempted to instantiate this class whenever `AKANTU_MODEL_COUPLERS` was enabled, leading to compilation failure.

**Modification**: Lines ~117-122: Added `#ifdef AKANTU_COHESIVE_ELEMENT` guard to ensure `CouplerSolidCohesiveContactOptions` is only used when the corresponding feature is compiled.

**Rationale**: This fix resolves a forward compatibility issue and follows proper conditional compilation practices for optional features.

### 2. Phase Field Model Enhancements (Research Contribution)

#### Motivation

Standard phase field fracture models typically employ a single energy functional that does not distinguish between tensile and compressive loading states. However, many quasi-brittle materials (concrete, ceramics, rocks, composites) exhibit markedly different behavior under tension versus compression: cracks readily propagate under tensile stress but remain closed under compression. This asymmetry is critical for accurately modeling fracture in structural materials and geomaterials.

The modifications implemented here introduce a volumetric-deviatoric energy split with tension-compression decomposition, enabling the phase field model to capture this physically realistic behavior. The compressive energy contribution (φ⁻) is tracked separately from the tensile contribution (φ⁺), allowing the degradation function to act primarily on tensile energy while preserving material stiffness under compression.

#### Implementation Details

**Core Data Structure Extensions** (`src/model/phase_field/phasefield.hh`, `phasefield.cc`)

Three new internal fields were introduced to the phase field model base class:
- `phi`: Phase field variable (damage parameter)
- `phi_minus`: Compressive component of elastic energy density (φ⁻)
- `phi_plus`: Tensile component of elastic energy density before history field tracking (φ⁺)

These fields are registered as internal variables in the finite element data structure and initialized to zero. Corresponding accessor methods (`getPhi()`, `getPhiMinus()`, `getPhiPlus()`) provide element-type-specific access to these fields for post-processing and visualization.

**Energy Split Formulation** (`energy_splits/volumetric_deviatoric_split_inline_impl.hh`)

The compressive energy contribution is computed at each quadrature point via the newly implemented `computePhiMinusOnQuad()` method:

φ⁻ = (1/2) κ ⟨tr(ε)⟩₋²

where:
- κ = λ + (2/3)μ is the bulk modulus
- tr(ε) is the volumetric strain (trace of strain tensor)
- ⟨·⟩₋ denotes the negative part: ⟨x⟩₋ = min(x, 0)
- λ, μ are the Lamé parameters

This formulation isolates the compressive volumetric energy, which remains undegraded by the phase field variable. A corresponding stub implementation was added to the `no_energy_split` class to maintain interface consistency when asymmetric splitting is disabled.

**Phase Field Model Integration** (`phasefields/phasefield_quadratic.cc`)

The quadratic phase field model was modified to compute and store both φ⁻ and φ⁺ during the energy computation phase. These values are used in the subsequent damage evolution calculation, where only the tensile energy drives crack growth.

**Python Interface Extensions** (`python/py_phasefield.cc`)

To facilitate post-processing and analysis of the asymmetric energy contributions, Python bindings were added for:
1. Field-specific accessors: `getPhi()`, `getPhiMinus()`, `getPhiPlus()`
2. Unified field access method: `getInternalReal(field_name, element_type)` supporting strain, phi, phi_minus, phi_plus, and damage fields

These bindings enable researchers to extract and visualize the energy decomposition directly from Python scripts, facilitating validation against analytical solutions and comparison with experimental data.

#### Theoretical Background

This implementation follows the volumetric-deviatoric split approach commonly employed in phase field models for quasi-brittle fracture [Amor et al., 2009; Miehe et al., 2010]. The key physical principle is that crack surfaces can only open under tension, not compression. By decomposing the elastic energy into tensile (φ⁺) and compressive (φ⁻) parts, and applying the degradation function g(d) only to φ⁺, the model prevents unphysical crack interpenetration under compressive loading while allowing crack propagation under tensile stress.

This approach has proven essential for modeling:
- Concrete structures under mixed-mode loading
- Hydraulic fracturing in rock formations
- Ceramic failure in thermal shock scenarios
- Composite delamination under bending

### 3. Documentation and Reproducibility

To ensure reproducibility and facilitate future development, comprehensive documentation has been prepared:

**`BUILD_MACOS_README.md`** - Detailed build guide including:
- Complete system specifications and dependency versions
- Step-by-step CMake configuration with all required flags
- Python bindings installation and troubleshooting
- Common build errors and their solutions
- Verification procedures

**`CREATE_BACKUP.md`** - Preservation strategies covering:
- Source code archival methods
- Binary distribution approaches
- Conda packaging for environment reproducibility
- Containerization for cross-platform deployment

**`README.md`** - Main repository documentation:
- Fork attribution and original source information
- Quick-start build instructions
- Summary of modifications
- Links to detailed documentation

These documents follow best practices for computational research reproducibility, enabling other researchers to rebuild this exact configuration on similar systems.

## Summary of Modifications

**Platform Compatibility**: 2 files modified (19 lines added)
- `cmake/Modules/FindMumps.cmake`: MUMPS detection on macOS
- `src/model/model.hh`: Conditional compilation fix

**Research Features**: 8 files modified (132 lines added, 1 removed)
- Phase field model base classes (2 files)
- Energy split implementations (4 files)
- Quadratic phase field model (1 file)
- Python bindings (1 file)

**Documentation**: 4 files created/modified
- BUILD_MACOS_README.md, CREATE_BACKUP.md, README.md, ATTRIBUTION.md

**Total Impact**: 13 files modified, 815 lines added, 1 line removed

## License and Copyright

This derived work maintains the original GNU Lesser General Public License v3.0 (LGPLv3) from Akantu. All modifications are contributed under the same license terms.

### Original Copyright
```
Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)

Akantu is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License along
with Akantu. If not, see <http://www.gnu.org/licenses/>.
```

## References

Key references for the phase field implementation approach:

1. Amor, H., Marigo, J. J., & Maurini, C. (2009). Regularized formulation of the variational brittle fracture with unilateral contact: Numerical experiments. *Journal of the Mechanics and Physics of Solids*, 57(8), 1209-1229.

2. Miehe, C., Hofacker, M., & Welschinger, F. (2010). A phase field model for rate-independent crack propagation: Robust algorithmic implementation based on operator splits. *Computer Methods in Applied Mechanics and Engineering*, 199(45-48), 2765-2778.

## Contact and Support

For questions regarding:
- **macOS build issues**: Open an issue on this repository
- **Phase field implementation**: Open an issue on this repository
- **General Akantu questions**: Refer to the original repository at https://gitlab.com/akantu/akantu
- **Original Akantu development**: Contact the LSMS laboratory at EPFL

---

*This work was developed as part of computational fracture mechanics research, building upon the excellent foundation provided by the Akantu development team at EPFL.*
