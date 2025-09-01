close all
addpath('./','./utils')

cf.C = 1;
cf.problem = 5;
cf.Nx1 = 217;
cf.Nx2 = 217;
cf.alpha = 1.01;

cf.nonlinear = 4;
cf.omega = 0.1;
cf.offdiag = 1;
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