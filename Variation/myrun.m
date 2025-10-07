close all
addpath('./','./utils')

cf.problem = 8;

cf.Nx1 = 41;
cf.Nx2 = 41;
cf.N = cf.Nx1*cf.Nx2;
cf.alpha = 1.005;

cf.sigma1 = 1;
cf.sigma2 = 1;

cf.nonlinear = 7;
cf.fixed_bc = 0;
cf.omega = 1;
cf.offdiag = 1;
cf.max_iter = 500;
cf.tolerance = 1e-6;
cf.smooth = 0;

cf.append_trail = 0;
cf.new_airfoil = 0;

cf.animation = 1;
cf.iter_per_frame = 1;
cf.plot_res = 1;
cf.pause_time = 0;
cf.make_gif = 0;
cf.gif_name = 'failed_example.gif';
cf.title_name = 'Alternative';
cf.save_output = 0;

if cf.problem == 1 || cf.problem == 2
    cf.C = 1;
end

if cf.nonlinear == 7
    cf.exact = 0;
    cf.repulse = 0;
end

if cf.new_airfoil == 1
    cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
    cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
end

[res_list,cost_list] = CurvedBv6(cf);