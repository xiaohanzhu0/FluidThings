# Recommended Layout

## Core API
- `+meshgen/` : new public package for mesh generation workflows.
  - `+meshgen/run_el.m` : entry point for EL formulation in physical coordinates.
  - `+meshgen/run_dual.m` : entry point for dual/harmonic-map formulation.
  - `+meshgen/defaults_el.m` : defaults + normalization for EL config (`cf`).
  - `+meshgen/defaults_dual.m` : defaults + normalization for dual config (`params`).

## Formulations
- `+meshgen/+formulations/`
  - `el_solve.m` : EL solver loop (refactored from `CurvedBv6.m`).
  - `dual_solve.m` : dual pipeline (refactored from `Harmonic/nonlinear_dual.m`).

## Problem & Metric Setup
- `+meshgen/+problems/`
  - `init_physical.m` : physical-grid initialization + boundary extraction.
  - `init_dual.m` : initialization for dual pipeline (grid + metric components).
- `+meshgen/+metrics/`
  - `for_problem.m` : metric function selection for EL formulation.

## Utilities
- `+meshgen/+util/`
  - `apply_defaults.m` : fill missing config fields with defaults.

## Legacy & Supporting Modules (unchanged)
- `+problems/` : problem definitions and metric helpers.
- `utils/` : linear system assembly, geometry utilities, I/O helpers.
- `Harmonic/` : harmonic map solver + dual assembly.
- `+analysis/` : mesh quality metrics and plotting.
- `test/regression/` : regression checks for refactor consistency.
