% set the boundary conditions from different minimization runs to use for
% LGM steady state run


% %OLD sets:
% load min1.mat
% load min2.mat
% load min3.mat
% load min4.mat
% load min5.mat
% load min6.mat
% load min7.mat
% load min8.mat
% load min9.mat  % this one combination of 4 (Darwin) and 6 (Hatherton)





% % -------------------------------------------------------------------------
% % If nothing is set in this file, defaults to load E and fs from
% % load_factors and load_factors2 in /TEST_DH
% % -------------------------------------------------------------------------


% % -------------------------------------------------------------------------

% % UNCOMMENT FOR STEP 1 and run snipet in turn to generate three files:
% % minE_values.mat; minfs_values.mat; minE_and_fs_values.mat
% % Sorry, this is an annoying step.

% %load minE_D.mat
% load minfs_D_2.mat
% %load minE_and_fs_D.mat
% x_P_min = x_P; x_w_min = x_w; x_e_min = x_e;
% E_P_min = E_P; E_w_min = E_w; E_e_min = E_e;
% fs_P_min = fs_P; fs_w_min = fs_w; fs_e_min = fs_e;
% 
% %load minE_H.mat
% load minfs_H.mat
% %load minE_and_fs_H.mat
% x_P2_min = x_P2; x_w2_min = x_w2; x_e2_min = x_e2;
% E_P2_min = E_P2; E_w2_min = E_w2; E_e2_min = E_e2;
% fs_P2_min = fs_P2; fs_w2_min = fs_w2; fs_e2_min = fs_e2;
% 
% save minfs_values.mat x_P_min x_w_min x_e_min E_P_min E_w_min E_e_min ...
%                      x_P2_min x_w2_min x_e2_min E_P2_min E_w2_min E_e2_min ...
%                      fs_P_min fs_w_min fs_e_min fs_P2_min fs_w2_min fs_e2_min

% % -------------------------------------------------------------------------


% % -------------------------------------------------------------------------

% % UNCOMMENT FOR STEP 2:

% % PICK ONE:
%load minE_values.mat
%load minfs_values.mat
load minE_and_fs_values.mat

% Interpolate to compare running on finer xgrid, whereas min search was
% done using "lower_resolution = 1"; interpolation can give negative (and
% then imaginary in S!) so just take absolute value as simple fix.
 E_P = abs(interp1(x_P_min, E_P_min, x_P, 'linear', 'extrap'));
 E_w = abs(interp1(x_w_min, E_w_min, x_w, 'linear', 'extrap'));
 E_e = abs(interp1(x_e_min, E_e_min, x_e, 'linear', 'extrap'));
 fs_P = abs(interp1(x_P_min, fs_P_min, x_P, 'linear', 'extrap'));
 fs_w = abs(interp1(x_w_min, fs_w_min, x_w, 'linear', 'extrap'));
 fs_e = abs(interp1(x_e_min, fs_e_min, x_e, 'linear', 'extrap'));

 E_P2 = abs(interp1(x_P2_min, E_P2_min, x_P2, 'linear', 'extrap'));
 E_w2 = abs(interp1(x_w2_min, E_w2_min, x_w2, 'linear', 'extrap'));
 E_e2 = abs(interp1(x_e2_min, E_e2_min, x_e2, 'linear', 'extrap'));
 fs_P2 = abs(interp1(x_P2_min, fs_P2_min, x_P2, 'linear', 'extrap'));
 fs_w2 = abs(interp1(x_w2_min, fs_w2_min, x_w2, 'linear', 'extrap'));
 fs_e2 = abs(interp1(x_e2_min, fs_e2_min, x_e2, 'linear', 'extrap'));


 
 