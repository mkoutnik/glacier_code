function residual = search_E ( factor_guess, ...
                               x_P, x_edges, x_w, x_e, dx_P, dx_w, dx_e, ...
                               B_P, B_w, B_e, W_P, W_w, W_e, ...
                               b_dot_edges, Q_0_in, S_0_in, ...
                               A_eff_edges_xt, BDM_surface, time, ...
                               measures_centerline_distance, ...
                               measures_flowspeed)
                                 


factor_w_guess = factor_guess(1:end-1);
factor_e_guess = factor_guess(2:end);


  [ S_P, h_P, ...
    dS_dx_P_xt,  ...
    dS_dx_edges_xt ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                 W_P, W_w, W_e, b_dot_edges(time,:), ...
                                 factor_w_guess, factor_e_guess, ...
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
%  flux_kin_P_xt(time,:) = calc_flux_kin( x_P, x_P, dx_P, W_P, ...
%                                         h_dot(time,:), b_dot_P(time,:), Q_0_in );
%                                                                                                  
% [ flux_w, flux_e ] = get_edge_values_quadratic ( flux_kin_P_xt(time,:), ...
%                                                  x_P, x_w, x_e, ...
%                                                  dx_P, dx_w, dx_e );
%    
%  flux_edges_kin_xt(time,:) = [ Q_0_in flux_e];   % Need to replace with Q_0_in, not extrapolated value.

                                                     
 flux_edges_dyn_xt(time,:) = calc_flux_dyn( x_edges, [h_w(time,1) h_e(time,:)], ...
                                            dS_dx_edges_xt(time,:), ...
                                            [factor_w_guess(1) factor_e_guess], ...
                                            [W_w(1) W_e], A_eff_edges_xt(time,:) );   
 flux_edges_dyn_xt(1,1)    = Q_0_in;                                    
          
 
 % Velocity
 % ========
 surf_vel_estimate = (5/4) * abs(flux_edges_dyn_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);

 BDM_velocity = interp1(measures_centerline_distance, measures_flowspeed, [x_w(1) x_e]);
 
 
 
 % Residual
 % ========
    std_dev = 1; 

 % % Surface profile:
 %  residual = sqrt( mean( ( (abs(BDM_surface) - abs(S_P))/std_dev ).^2 ) );           

 % Surface velocity:
   residual = sqrt( mean( ( (abs(BDM_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
   
   
   figure(1)
   plot(x_edges, BDM_velocity, 'k', 'linewidth', 2)
   hold on
   plot(x_edges, surf_vel_estimate, 'r')
   
   figure(2)
   plot(x_edges, [factor_w_guess(1) factor_e_guess],'b')
   hold on
   
   
   

