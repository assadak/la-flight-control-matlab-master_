clear all
close all
clc

addpath(fullfile('3rd_party', 'YAMLMatlab_0.4.3'))
addpath(genpath('src'))
 
approx = PolyApprox(2.5, 2, 2);
% approx.approx([1, 2], 0.5)
% cgl1 = approx.cgl_points(4, 0, 1);
% cheb_points = approx.cheb_points(3, 0, 1)

approx.approx_test()


%%%