# Changelog

## Version 1.0.0 (2026-04-08)

### Features
- **Generalized Qn Parameter**: Support for any Steinhardt order parameter (Q1-Q11+)
- **Flexible Pair Analysis**: Compute different order parameters for different pairs
- **Two Implementations**:
  - V0: Sequential version for debugging and single pairs
  - V3: Vectorized + parallelized version for production
- **Center-of-Mass Support**: Optional COM calculation for molecular systems
- **Adaptive Cutoff**: Automatic or manual neighbor distance cutoff
- **Detailed Logging**: Comprehensive neighbor information and statistics
- **Multi-pair Tracking**: Independent coupling per pair

### Architecture
- **Performance**: V3 achieves 6-12x speedup over V0 for multi-pair analysis
- **Vectorization**: NumPy cdist() and batch spherical harmonic computation
- **Parallelization**: Frame-level multiprocessing with worker pool
- **Memory Efficient**: Scales well from small to large systems

### Output Files
- `Qn.xvg` - Multi-column time series
- `Qn_NN.xvg` - Neighbor count tracking
- `Qn.log` - Detailed statistics and execution log
- `Qn_detail.txt` - Optional neighbor lists

### Documentation
- Comprehensive README with examples
- Technical MANUAL with algorithm details
- Quick start guide
- Example input files and use cases

### Known Limitations
- Requires GROMACS trajectory in .gro format
- GPU acceleration not yet implemented
- Bonding topology not used for neighbor detection

### Future Work
- GPU CUDA implementation
- Weighted harmonic computation
- Symmetry parameter (Wl) calculation
- Bond topology integration

## Version 0.9.0 (Development)

### Initial Development
- Q6 implementation with all core features
- Neighbor deduplication algorithm
- Molecule-level self-pair exclusion
- V0 and V3 versions with identical output verification
