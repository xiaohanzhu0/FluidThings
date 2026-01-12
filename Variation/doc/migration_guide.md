# Migration Guide

## Entry Points
- `myrun.m` -> `meshgen.run_el` (same `cf` fields). `myrun.m` still works as a legacy wrapper.
- `CurvedBv6.m` -> `meshgen.formulations.el_solve` (now returns `[x1, x2, res_list, cost_list, info]`).
- `Harmonic/nonlinear_dual.m` -> `meshgen.run_dual` or `meshgen.formulations.dual_solve`.

## Function Relocations
- EL solver loop: `CurvedBv6.m` -> `+meshgen/+formulations/el_solve.m`.
- Dual pipeline: `Harmonic/nonlinear_dual.m` -> `+meshgen/+formulations/dual_solve.m`.
- Metric selection: `utils/GetM.m` -> `+meshgen/+metrics/for_problem.m` (wrapper).
- Physical-grid initialization: inline logic in `CurvedBv6.m` -> `+meshgen/+problems/init_physical.m`.
- Dual initialization: `+problems/Initialization.m` + `utils/extract_boundary_points.m` -> `+meshgen/+problems/init_dual.m`.
- Default configurations: `+meshgen/defaults_el.m`, `+meshgen/defaults_dual.m`.

## Notes
- Legacy entry points remain available and call the new implementations.
- Plotting is controlled by flags: `cf.animation`, `cf.plot_res`, `cf.make_gif` (EL) and `params.doPlot`, `params.plotInitial`, `params.showHarmonicPlot`, `params.showHarmonicRes` (dual).
- Metric gradation/smoothing behavior is unchanged; see `+problems/InitProb5.m` and `utils/grade_structured_mesh.m`.
