% run_min_search_E

E_vec         = [ 0.5:0.5:20 ];
N_E           = length( E_vec );

length_run    = length(x_P);
save_E        = [];
save_mismatch = NaN * ones(length_run, N_E);

factor_P_orig = factor_P;


for xpos = length_run:-1:1 % 1:length_run   % loop over all x positions

xpos                       % display to screen

RMS_mismatch = NaN * ones( 1, N_E );


%  Loop over E_x
 for i_E = 1:N_E

      E_use = E_vec(i_E);

      factor_P(xpos) = E_use;
      
    % recalculate values on the edges.
      [factor_w, factor_e ] = get_edge_values_quadratic ...
                               ( factor_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

  [ S_P, h_P, ...
    dS_dx_P_xt,  ...
    dS_dx_edges_xt ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                 W_P, W_w, W_e, b_dot_edges(time,:), ...
                                 factor_w, factor_e, ...
                                 Q_0_in, S_0_in, ...
                                 A_eff_edges_xt(time,1:end-1), ...
                                 A_eff_edges_xt(time,2:end) );
                                     
% edge values of ice thickness
% ============================
  [ h_w(time,:), ...
    h_e(time,:) ] = get_edge_values_quadratic( h_P(time,:), x_P, x_w, x_e, ...
                                               dx_P, dx_w, dx_e );
              
% Find the flux
% ==============
 flux_edges_dyn_xt(time,:) = calc_flux_dyn( x_edges, [h_w(time,1) h_e(time,:)], ...
                                            dS_dx_edges_xt(time,:), ...
                                            [factor_w(1) factor_e], ...
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
        
   RMS_mismatch(i_E) = residual;
   
                                   
 end  % loop over E values
 
 
   index     = find( RMS_mismatch == min(min(RMS_mismatch) ) );
   RMS_best  = RMS_mismatch(index)
   best_E    = E_vec(index)
   save_E    = [ save_E best_E ];    
   save_mismatch(xpos, :) = RMS_mismatch;

% reset E:
  factor_P(xpos) = best_E;
 
end  % loop over xpositions


 % recalculate values on the edges.
      [factor_w, factor_e ] = get_edge_values_quadratic ...
                               ( factor_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

 [ S_P(time,:), h_P(time,:), ...
    dS_dx_P_xt(time,:),  ...
    dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                         W_P, W_w, W_e, b_dot_edges(time,:), ...
                                         factor_w, factor_e, ...
                                         Q_0_in, S_0_in, ...
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
                                            [factor_w(1) factor_e], ...
                                            [W_w(1) W_e], A_eff_edges_xt(time,:) );   
 flux_edges_dyn_xt(1,1)    = Q_0_in; 

 
 average_vel_estimate = abs(flux_edges_dyn_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);

                           
                           

figure(10)
plot(x_P, factor_P_orig,'k')
hold on
plot(x_P, factor_P, 'r--')
title('Enhancement factor, E')


% -------------------------------------------------------------------------