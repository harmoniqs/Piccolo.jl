# Release Notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.0] - Upcoming

### Changed

- Complete documentation overhaul with concept-driven organization
- Migrated from PiccoloDocsTemplate to standard Documenter.jl

## [v0.3.1] - 2024-10-17

### Fixed

- Fixed and added tests to `RydbergChainSystem`

## [v0.3.0] - 2024-10-10

### Added

- `PiccoloOptions` to handle custom problem settings

### Changed

- Refactored trajectory initialization functions
- Improved documentation
- Typo fixes

## [v0.2.0] - 2024-02-22

### Added

- `EmbeddedOperator` to handle subspace gate optimization and leakage suppression
- Plotting methods for unitary populations

### Changed

- New quantum systems interface
  - Transmon system template
- Restructured the code base for easier quantum system and problem template development

### Removed

- Stale examples

### Fixed

- Robustness improvements objective test fixes
