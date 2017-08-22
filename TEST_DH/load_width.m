function [ W_P, W_w, W_e ] = load_width( x_P, x_w, x_e, dx_P, dx_w, dx_e )


%--------------------------------------------------------------
%
% Michelle Koutnik (mkoutnik@ess.washington.edu)
%
%  set up steady width values W_P(x) at positions x_P
%  where surface S(x,t) is defined
%  interpolate bed to interval midpoints where fluxes are
%  calculated
%
%---------------------------------------------------------------

global DIRECTORY_data

addpath(DIRECTORY_data)

               
% % ------------- 
% % Uniform width
% % -------------
% 
%     width = 1;     
%     W_P = width * ones( size(x_P) );

    
% --------------- 
% Measured width
% ---------------

% % From QGIS -- width
load DH_accum_width_velocity.mat
W_P = interp1(Darwin_width_x, Darwin_width_values, x_P);
W_P = smooth(W_P./W_P(1),2)';


% % Values from min search
% % ======================
% disp('Using best min search width!')
% load best_min_search_width.mat
% W_P = interp1(x_P_best, W_P_best, x_P, 'linear', 'extrap');

  

% need to extrapolate to control-volume edges and interpolate for
% the interior points

  [W_w, W_e ] = get_edge_values_quadratic ...
                               ( W_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );


