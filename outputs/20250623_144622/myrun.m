addpath('./','./utils')

cf.C = 1;
cf.nonlinear = 4;
cf.omega = 0.5;
cf.problem = 1;
cf.Nx1 = 40;
cf.Nx2 = 40;
cf.offdiag = 0;
cf.max_iter = 500;
cf.tolerance = 1e-6;

cf.animation = 1;
cf.plot_res = 1;
cf.pause_time = 0;

cf.make_gif = 0;
cf.gif_name = 'example13.gif';
cf.title_name = 'Sampled airfoil, sampled metric';
cf.save_output = 1;

[res_list,cost_list] = CurvedBv6(cf);