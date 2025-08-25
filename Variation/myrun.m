close all
addpath('./','./utils')

cf.C = 2;
cf.problem = 3;
cf.Nx1 = 217*2+1;
cf.Nx2 = 71*1;
cf.Nx1 = 40;
cf.Nx2 = 40;
cf.N = cf.Nx1*cf.Nx2;
cf.alpha = 1.005;

cf.nonlinear = 7;
cf.fixed_bc = 0;
cf.exact = 0;
cf.omega = 1;
cf.offdiag = 1;
cf.repulse = 0;
cf.max_iter = 500;
cf.tolerance = 1e-9;
cf.smooth = 0;

cf.append_trail = 0;
cf.new_airfoil = 1;
cf.animation = 1;
cf.iter_per_frame = 1;
cf.plot_res = 1;
cf.pause_time = 0;

cf.make_gif = 0;
cf.gif_name = 'failed_example.gif';
cf.title_name = 'Alternative';
cf.save_output = 0;

[res_list,cost_list] = CurvedBv6(cf);