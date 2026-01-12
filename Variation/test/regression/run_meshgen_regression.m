function results = run_meshgen_regression()
%RUN_MESHGEN_REGRESSION Compare new entry points against core solvers.
    tol = 1e-10;
    el_grid = 21;
    dual_grid = 21;
    dual_sample = 31;

    results = struct('el', [], 'dual', []);

    fprintf('Running EL formulation regression...\n');
    for problemId = 1:5
        cf = meshgen.defaults_el(struct(...
            'problem', problemId, ...
            'Nx1', el_grid, ...
            'Nx2', el_grid, ...
            'animation', 0, ...
            'plot_res', 0, ...
            'make_gif', 0, ...
            'save_output', 0, ...
            'pause_time', 0 ...
        ));

        if problemId == 5 && ~airfoil_files_available(cf.metric_datapath, cf.airfoil_datapath)
            fprintf('  EL problem %d skipped (airfoil data missing).\n', problemId);
            results.el(problemId).skipped = true; %#ok<*AGROW>
            continue;
        end

        [x1_core, x2_core] = meshgen.formulations.el_solve(cf);
        [x1_entry, x2_entry] = meshgen.run_el(cf);

        [dx1, dx2] = compare_meshes(x1_core, x2_core, x1_entry, x2_entry, tol);
        results.el(problemId).max_dx1 = dx1;
        results.el(problemId).max_dx2 = dx2;
        results.el(problemId).skipped = false;
        fprintf('  EL problem %d: max|dx1|=%.3e, max|dx2|=%.3e\n', problemId, dx1, dx2);
    end

    fprintf('Running dual formulation regression...\n');
    for problemId = 1:5
        params = meshgen.defaults_dual(struct(...
            'problemId', problemId, ...
            'Nx1', dual_grid, ...
            'Nx2', dual_grid, ...
            'Nt1', dual_grid, ...
            'Nt2', dual_grid, ...
            'sampleNs', dual_sample, ...
            'doPlot', 0, ...
            'plotInitial', 0, ...
            'showHarmonicPlot', false, ...
            'showHarmonicRes', false ...
        ));

        if problemId == 5
            params_el = meshgen.defaults_el(struct('problem', 5));
            if ~airfoil_files_available(params_el.metric_datapath, params_el.airfoil_datapath)
                fprintf('  Dual problem %d skipped (airfoil data missing).\n', problemId);
                results.dual(problemId).skipped = true;
                continue;
            end
        end

        [x1_core, x2_core] = meshgen.formulations.dual_solve(params);
        [x1_entry, x2_entry] = meshgen.run_dual(params);

        [dx1, dx2] = compare_meshes(x1_core, x2_core, x1_entry, x2_entry, tol);
        results.dual(problemId).max_dx1 = dx1;
        results.dual(problemId).max_dx2 = dx2;
        results.dual(problemId).skipped = false;
        fprintf('  Dual problem %d: max|dx1|=%.3e, max|dx2|=%.3e\n', problemId, dx1, dx2);
    end
end

function [dx1, dx2] = compare_meshes(x1_ref, x2_ref, x1_test, x2_test, tol)
    dx1 = max(abs(x1_ref(:) - x1_test(:)));
    dx2 = max(abs(x2_ref(:) - x2_test(:)));
    if dx1 > tol || dx2 > tol
        error('Regression check failed: max|dx1|=%.3e, max|dx2|=%.3e', dx1, dx2);
    end
end

function tf = airfoil_files_available(metric_path, airfoil_path)
    metric_path = expand_user(metric_path);
    airfoil_path = expand_user(airfoil_path);
    tf = exist(metric_path, 'file') == 2 && exist(airfoil_path, 'file') == 2;
end

function path_out = expand_user(path_in)
    if startsWith(path_in, '~/')
        home_dir = getenv('HOME');
        path_out = fullfile(home_dir, path_in(3:end));
    elseif startsWith(path_in, '~')
        home_dir = getenv('HOME');
        path_out = fullfile(home_dir, path_in(2:end));
    else
        path_out = path_in;
    end
end
