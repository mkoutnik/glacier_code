% calculate_surface.m    


% ------------------------------------------------------------------------

% script with all calls to .m files 
%
% ------------------------------------------------------------------------


% Assume the temperature field, then calculate the ice-surface evolution
% and the velocity field over time

global min_search_E min_search_fs min_search_bed
global figures_on
global DIRECTORY_data

addpath(DIRECTORY_data)





% NEED TO HAVE FILE ONHAND TO COMPARE TO! First run with all min_search = 0
% and steady_state = 1
if (figures_on == 1)
    
load values_nochanges.mat    
    
figure(1)   % SURFACE
plot(x_P, B_P,'r')
hold on
plot(x_P, S_modern)
plot(x_P_nochanges, S_P_nochanges(1,:),'k')
plot(x_P, S_P(1,:),'c--')
legend('Bed', 'Obs surface', 'Calculated surface no changes', 'Calculated surface', 'location', 'northwest')

% figure(2)   % FLUX
% plot(x_edges, flux_edges_kin_xt(1,:),'b')
% hold on
% plot(x_edges, flux_edges_dyn_xt(1,:),'r--')
% %plot(x_edges, flux_edges_kin_xt(1,:) - flux_edges_dyn_xt(1,:))
        
figure(3)   % SURFACE VELOCITY
plot(x_edges_nochanges, (5/4)* average_vel_estimate_nochanges,'k')
hold on
plot(x_edges, (5/4)* average_vel_estimate,'c--')
plot(measures_centerline_distance, measures_flowspeed,'m', 'linewidth', 2)
legend('Calculated no changes','Calculated',  'MEaSUREs values', 'location', 'best')
xlim([0 x_P(end-1)])
title('Surface velocity')

% figure(4)   % THICKNESS
% hold on
% plot(x_P, S_modern-B_P)
% plot(x_P_nochanges, S_P_nochanges(1,:)-B_P,'k')
% plot(x_P, S_P(1,:)-B_P,'c--')
% legend('Obs surface', 'Calculated surface no changes', 'Calculated surface', 'location', 'northwest')
% title('Ice thickness')

figure(5)   % THICKNESS DIFFERENCE
plot(x_P, (S_modern-B_P) - (S_P(1,:)-B_P))
title('Thickness difference (h_obs - h_calc)')

end




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    

% ----------------        
if (time == 1)    
% ----------------    
    
% Guess initial surface...    
    
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
                                        
 
 average_vel_estimate = abs(flux_edges_dyn_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
                                        
                                        
   
    

% assign steady-state values
% ==========================
   Q_out_L_SS     = flux_edges_dyn_xt(time,1);
   Q_out_R_SS     = flux_edges_dyn_xt(time,end);
 % these values are used in volume_perturbation.m for boundary condition                      
  
   Q_out_L(time)  = Q_out_L_SS;
   Q_out_R(time)  = Q_out_R_SS;
 
  
  
   
% ----------------    
elseif (time > 1)  % this comes in if steady_state_only = 0
% ----------------    
    
    

% initialize values used in the surface solver for the future timestep to 
% those at the previous timestep for the first iteration.
% =======================================================
   A_eff_edges_xt(time,:) = A_eff_edges_xt(time-1,:);
   u_bar_edges_xt(time,:) = u_bar_edges_xt(time-1,:);
   S_P(time,:)            = S_P(time-1,:);
   h_P(time,:)            = h_P(time-1,:);
   S_dot(time,:)          = S_dot(time-1,:);
   T_field_xzt(time,:,:)  = T_field_xzt(time-1,:,:);
   
    

% --------------------------------
% TRANSIENT SURFACE CALCULATION 
% --------------------------------

DIRECTORY_surf = 'MFILES_surface_main';

addpath( DIRECTORY_surf )   

% % always send the known values at time-1 and estimates of values at time
% % (which in the first iteration are the values at time-1)
                             
   
% This is clunky... steady-state solver and first timestep of transient solver
% don't give the same... does difference come in because of spatial step? What gives?
 if (time == 2)    
  hold_SS_check = 9999;
 end
 while (hold_SS_check > 0.01)  % Surface profile not holding SS to 0.01 m
     
      disp(' checking holds steady state...')
      [ h_P_check, h_w_check, h_e_check, ...
        S_P_check, S_w_check, S_e_check, ...
        flux_edges_dyn_xt_check ] = surface_main( h_P(time-1,:), S_P(time-1,:), ...
                                                  x_P, x_w, x_e, ...
                                                  dx_P, dx_w, dx_e, ...
                                                  B_P, B_w, B_e, W_P, W_w, W_e,  ...
                                                  E_P, E_w, E_e, ...
                                                  fs_P, fs_w, fs_e, ...
                                                  A_eff_edges_xt(time-1:time,:), ...
                                                  b_dot_P(time-1:time,:), ...
                                                  b_dot_edges_SS, ...
                                                  b_dot_edges(time-1:time,:), ...
                                                  b_dot_P, b_dot_edges, ...
                                                  Q_out_L_SS, Q_out_R_SS, ...
                                                  S_at_GL(time-1:time), ...
                                                  Q_ext_P(time-1:time,1), ...
                                                  Q_ext_P(time-1:time,2), ...  
                                                  dt_P(time), t_P, time );
     
   hold_SS_check = max (S_P(time,:) - S_P_check);                                   
                                          
   h_P(time,:) = h_P_check;
   h_w(time,:) = h_w_check;
   h_e(time,:) = h_e_check;
   S_P(time,:) = S_P_check;
   S_w(time,:) = S_w_check;
   S_e(time,:) = S_e_check;
   flux_edges_dyn_xt(time,:) = flux_edges_dyn_xt_check;
   
 end
 
 
  [ h_P(time,:), h_w(time,:), h_e(time,:), ...
   S_P(time,:), S_w(time,:), S_e(time,:), ...
   flux_edges_dyn_xt(time,:)      ] = surface_main( h_P(time-1,:), S_P(time-1,:), ...
                                              x_P, x_w, x_e, ...
                                              dx_P, dx_w, dx_e, ...
                                              B_P, B_w, B_e, W_P, W_w, W_e,  ...
                                              E_P, E_w, E_e, ...
                                              fs_P, fs_w, fs_e, ...
                                              A_eff_edges_xt(time-1:time,:), ...
                                              b_dot_P(time-1:time,:), ...
                                              b_dot_edges_SS, ...
                                              b_dot_edges(time-1:time,:), ...
                                              b_dot_P, b_dot_edges, ...
                                              Q_out_L_SS, Q_out_R_SS, ...
                                              S_at_GL(time-1:time), ...
                                              Q_ext_P(time-1:time,1), ...
                                              Q_ext_P(time-1:time,2), ...  
                                              dt_P(time), t_P, time );
            
                                          
% fill matrices to save values
% ============================
     S_P(time,:)                = h_P(time,:) + B_P; 
    
     h_dot(time,:)              = (h_P(time,:) - h_P(time-1,:)) / dt_P(time);
     
     S_dot(time,:)              = (S_P(time,:) - S_P(time-1,:)) / dt_P(time);
     
     [ dS_dx_w, dS_dx_e ]       = get_gradient_values( S_P(time,:), x_P, dx_P );
     dS_dx_edges_xt(time,:)     = [dS_dx_w(1) dS_dx_e];
     
     dS_dx_P_xt(time,:)         = qinterp1(x_edges, dS_dx_edges_xt(time,:), x_P)';
    
     
% FLUX
% ====
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

 
 average_vel_estimate = abs(flux_edges_dyn_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
                                        
                                        

rmpath( DIRECTORY_surf )    % remove path of MFILES_surface_main    



end   % if statement on time



z_grid = ( (z_hat_grid' .* repmat(h_P(time,:), N_z, 1)) + repmat(B_P, N_z, 1) )';


% ____________________________________________________   end loop on time     



% Save file with no parameter changes to compare against in figures:
if ((min_search_bed == 0) && (min_search_E == 0) && (min_search_fs == 0))
    
x_P_nochanges                  = x_P;
S_P_nochanges                  = S_P;
x_edges_nochanges              = x_edges;
average_vel_estimate_nochanges = average_vel_estimate;

save values_nochanges.mat x_P_nochanges S_P_nochanges x_edges_nochanges average_vel_estimate_nochanges

end




