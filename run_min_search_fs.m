% run_min_search_fs



fs_vec     = [ 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 ]; % 1e-8 1e-7 1e-6 1e-5 ];
N_fs       = length( fs_vec );

length_run    = length(x_P);
save_fs       = [];
save_mismatch = NaN * ones(length_run, N_fs);

fs_P_orig = fs_P;


for xpos = length_run:-1:1 % 1:length_run   % loop over all x positions

 disp(['Evaluating fs for x-position ',int2str(xpos), ' of total=',int2str(length_run)]);

RMS_mismatch = NaN * ones( 1, N_fs );


%  Loop over fs values
 for i_fs = 1:N_fs
 
      fs_P   = fs_P_orig;
      fs_use = fs_vec(i_fs);

      fs_P(xpos) = fs_use;
      
    % recalculate values on the edges.
      [fs_w, fs_e ] = get_edge_values_quadratic ...
                               ( fs_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
                           
  [ S_P(time,:), h_P(time,:), ...  
    dS_dx_P_xt(time,:),  ...
    dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                         W_P, W_w, W_e, b_dot_edges(time,:), ...
                                         E_w, E_e, fs_w, fs_e, ...
                                         Q_0_in, S_at_GL(time), ...
                                         A_eff_edges_xt(time,1:end-1), ...
                                         A_eff_edges_xt(time,2:end) );
                                     
                         
% edge values of ice thickness
% ============================
  [ h_w(time,:), ...
    h_e(time,:) ] = get_edge_values_quadratic( h_P(time,:), x_P, x_w, x_e, ...
                                               dx_P, dx_w, dx_e );
              
% Find the flux
% ==============
 flux_kin_P_xt(time,:) = calc_flux_kin( x_P, x_P, dx_P, W_P, ...
                                        h_dot(time,:), b_dot_P(time,:), Q_0_in );
                                                                                             
[ flux_w, flux_e ] = get_edge_values_quadratic ( flux_kin_P_xt(time,:), ...
                                                 x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );
   
 flux_edges_kin_xt(time,:) = [ Q_0_in flux_e];   % Need to replace with Q_0_in, not extrapolated value.
                                                    
 flux_edges_dyn_xt(time,:) = calc_flux_dyn( x_edges, [h_w(time,1) h_e(time,:)], ...
                                            dS_dx_edges_xt(time,:), ...
                                            [E_w(1) E_e], [fs_w(1) fs_e], ...
                                            [W_w(1) W_e], A_eff_edges_xt(time,:) );   
 flux_edges_dyn_xt(1,1)    = Q_0_in; 
                                        
           
 
 % Velocity
 % ========
 surf_vel_estimate = (5/4) * abs(flux_edges_dyn_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
 
 BDM_velocity      = interp1(measures_centerline_distance, measures_flowspeed, [x_w(1) x_e]);
 
    
 
 
 % Residual
 % ========
    std_dev = 1; 

 % % Surface profile:
 %  residual = sqrt( mean( ( (abs(BDM_surface) - abs(S_P))/std_dev ).^2 ) );           

 % Surface velocity:
   residual = sqrt( mean( ( (abs(BDM_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
        
   RMS_mismatch(i_fs) = residual;
   
                                   
 end  % loop over fs values
 
 
   index      = find( RMS_mismatch == min(min(RMS_mismatch) ) );
   RMS_best   = RMS_mismatch(index);
   best_fs    = fs_vec(index);
   save_fs    = [ save_fs best_fs ];    
   save_mismatch(xpos, :) = RMS_mismatch;

% % reset fs_P:
%  fs_P(xpos) = best_fs;
 
end  % loop over xpositions



 fs_P = save_fs;
      
    % recalculate values on the edges.
      [ fs_w, fs_e ] = get_edge_values_quadratic ...
                               ( fs_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );




% -------------------------------------------------------------------------