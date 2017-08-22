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





% % NEED TO HAVE FILE ONHAND TO COMPARE TO! First run with all min_search = 0
% % and steady_state = 1
% if (figures_on == 1)
%     
% load values_nochanges.mat    
%     
% figure(1)   % SURFACE
% plot(x_P, B_P,'r')
% hold on
% plot(x_P, S_modern)
% plot(x_P_nochanges, S_P_nochanges(1,:),'k')
% plot(x_P, S_P(1,:),'c--')
% legend('Bed', 'Obs surface', 'Calculated surface no changes', 'Calculated surface', 'location', 'northwest')
% 
% % figure(2)   % FLUX
% % plot(x_edges, flux_edges_kin_xt(1,:),'b')
% % hold on
% % plot(x_edges, flux_edges_dyn_xt(1,:),'r--')
% % %plot(x_edges, flux_edges_kin_xt(1,:) - flux_edges_dyn_xt(1,:))
%         
% figure(3)   % SURFACE VELOCITY
% plot(x_edges_nochanges, (5/4)* average_vel_estimate_nochanges,'k')
% hold on
% plot(x_edges, (5/4)* average_vel_estimate,'c--')
% plot(measures_centerline_distance, measures_flowspeed,'m', 'linewidth', 2)
% legend('Calculated no changes','Calculated',  'MEaSUREs values', 'location', 'best')
% xlim([0 x_P(end-1)])
% title('Surface velocity')
% 
% % figure(4)   % THICKNESS
% % hold on
% % plot(x_P, S_modern-B_P)
% % plot(x_P_nochanges, S_P_nochanges(1,:)-B_P,'k')
% % plot(x_P, S_P(1,:)-B_P,'c--')
% % legend('Obs surface', 'Calculated surface no changes', 'Calculated surface', 'location', 'northwest')
% % title('Ice thickness')
% 
% figure(5)   % THICKNESS DIFFERENCE
% plot(x_P, (S_modern-B_P) - (S_P(1,:)-B_P))
% title('Thickness difference (h_obs - h_calc)')
% 
% end




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    

% ----------------        
if (time2 == 1)    
% ----------------    
    
% Guess initial surface...    
    
   [ S_P2(time2,:), h_P2(time,:), ...  
    dS_dx_P_xt2(time2,:),  ...
    dS_dx_edges_xt2(time2,:) ] = calc_h_0( x_P2, x_w2, x_e2, dx_P2, B_P2, B_w2, B_e2, ...
                                         W_P2, W_w2, W_e2, b_dot_edges2(time2,:), ...
                                         E_w2, E_e2, fs_w2, fs_e2, ...
                                         Q_0_in2, S_at_GL2(time2), ...
                                         A_eff_edges_xt2(time2,1:end-1), ...
                                         A_eff_edges_xt2(time2,2:end), ...
                                         flux_add_w2, flux_add_e2, ...
                                         deformation_only2, deformation_plus_sliding2, sliding_only2);
                                     
                         
% edge values of ice thickness
% ============================
  [ h_w2(time2,:), ...
    h_e2(time2,:) ] = get_edge_values_quadratic( h_P2(time2,:), x_P2, x_w2, x_e2, ...
                                                 dx_P2, dx_w2, dx_e2 );
              
% Find the flux
% ==============
 flux_kin_P_xt2(time2,:) = calc_flux_kin( x_P2, x_P2, dx_P2, W_P2, ...
                                        h_dot2(time2,:), b_dot_P2(time2,:), Q_0_in2 );
                                                                                             
[ flux_w2, flux_e2 ] = get_edge_values_quadratic ( flux_kin_P_xt2(time2,:), ...
                                                 x_P2, x_w2, x_e2, ...
                                                 dx_P2, dx_w2, dx_e2 );
   
 flux_edges_kin_xt2(time2,:) = [ Q_0_in2 flux_e2];   % Need to replace with Q_0_in, not extrapolated value.
                                               
 
%  flux_edges_dyn_xt2(time2,:) = calc_flux_dyn( x_edges2, [h_w2(time2,1) h_e2(time2,:)], ...
%                                             -dS_dx_edges_xt2(time2,:), ...
%                                             [E_w2(1) E_e2], [fs_w2(1) fs_e2], ...
%                                             [W_w2(1) W_e2], A_eff_edges_xt2(time2,:) );   
% 
                                                     
 flux_edges_dyn_xt2(time2,:) = calc_flux_dyn( x_edges2, [h_w2(time2,1) h_e2(time2,:)], ...
                                            dS_dx_edges_xt2(time2,:), ...
                                            [E_w2(1) E_e2], [fs_w2(1) fs_e2], ...
                                            [W_w2(1) W_e2], A_eff_edges_xt2(time2,:), ...
                                            deformation_only2, deformation_plus_sliding2, sliding_only2);   
 
 % % CHECK THIS.                                       
 flux_edges_dyn_xt2(1,1)    = Q_0_in; 
                 
 flux_edges_dyn_xt2(1,:) = flux_edges_dyn_xt2(1,:) + [flux_add_w2(1) flux_add_e2];
 
 
% average_vel_estimate2 = abs(flux_edges_dyn_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
                                        
 average_vel_estimate2 = abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
                                        

%  plot(flux_edges_dyn_xt2(1,:),'c')
%  hold on
%  plot(flux_edges_kin_xt2(1,:),'r--')

    

% assign steady-state values
% ==========================
   Q_out_L_SS2     = flux_edges_kin_xt2(time2,1);
   Q_out_R_SS2     = flux_edges_kin_xt2(time2,end);
 % these values are used in volume_perturbation.m for boundary condition                      
  
   Q_out_L2(time2)  = Q_out_L_SS2;
   Q_out_R2(time2)  = Q_out_R_SS2;
 
  
  
   
% ----------------    
elseif (time > 1)  % this comes in if steady_state_only = 0
% ----------------    
    
    

% initialize values used in the surface solver for the future timestep to 
% those at the previous timestep for the first iteration.
% =======================================================
   A_eff_edges_xt2(time2,:) = A_eff_edges_xt2(time2-1,:);
   u_bar_edges_xt2(time2,:) = u_bar_edges_xt2(time2-1,:);
   S_P2(time2,:)            = S_P2(time2-1,:);
   h_P2(time2,:)            = h_P2(time2-1,:);
   S_dot2(time,:)           = S_dot2(time2-1,:);
   T_field_xzt2(time2,:,:)  = T_field_xzt2(time2-1,:,:);
   
    

% --------------------------------
% TRANSIENT SURFACE CALCULATION 
% --------------------------------

DIRECTORY_surf = 'MFILES_surface_main';

addpath( DIRECTORY_surf )   

% % always send the known values at time-1 and estimates of values at time
% % (which in the first iteration are the values at time-1)
                             
   
% % This is clunky... steady-state solver and first timestep of transient solver
% % don't give the same... does difference come in because of spatial step? 
% % Because dynamic and kinematic flux off depending on resolution? What gives?
%  if (time2 == 2)    
%   hold_SS_check = 9999;
%  end
%  while (hold_SS_check > 0.001)  % Surface profile not holding SS to 0.01 m
%      
%       disp(' checking holds steady state...')
%       [ h_P_check2, h_w_check2, h_e_check2, ...
%         S_P_check2, S_w_check2, S_e_check2, ...
%         flux_edges_dyn_xt_check2 ] = surface_main( h_P2(time2-1,:), S_P2(time2-1,:), ...
%                                                   x_P2, x_w2, x_e2, ...
%                                                   dx_P2, dx_w2, dx_e2, ...
%                                                   B_P2, B_w2, B_e2, W_P2, W_w2, W_e2,  ...
%                                                   E_P2, E_w2, E_e2, ...
%                                                   fs_P2, fs_w2, fs_e2, ...
%                                                   A_eff_edges_xt2(time2-1:time2,:), ...
%                                                   b_dot_P2(time-1:time,:), ...
%                                                   b_dot_edges_SS2, ...
%                                                   b_dot_edges2(time2-1:time2,:), ...
%                                                   b_dot_P2, b_dot_edges2, ...
%                                                   Q_out_L_SS2, Q_out_R_SS2, ...
%                                                   S_at_GL2(time2-1:time2), ...
%                                                   Q_ext_P2(time2-1:time2,1), ...
%                                                   Q_ext_P2(time2-1:time2,2), ...  
%                                                   dt_P2(time2), t_P2, time2, ...
%                                                   deformation_only2, deformation_plus_sliding2, sliding_only2);
%      
%    hold_SS_check = max (S_P2(time2,:) - S_P_check2)                                
%                                           
%    h_P2(time2,:) = h_P_check2;
%    h_w2(time2,:) = h_w_check2;
%    h_e2(time2,:) = h_e_check2;
%    S_P2(time2,:) = S_P_check2;
%    S_w2(time2,:) = S_w_check2;
%    S_e2(time2,:) = S_e_check2;
%    flux_edges_dyn_xt2(time2,:) = flux_edges_dyn_xt_check2;
%    
%  end
 
 
  [ h_P2(time2,:), h_w2(time2,:), h_e2(time2,:), ...
   S_P2(time2,:), S_w2(time2,:), S_e2(time2,:), ...
   flux_edges_dyn_xt2(time2,:)      ] = surface_main( h_P2(time2-1,:), S_P2(time2-1,:), ...
                                              x_P2, x_w2, x_e2, ...
                                              dx_P2, dx_w2, dx_e2, ...
                                              B_P2, B_w2, B_e2, W_P2, W_w2, W_e2,  ...
                                              E_P2, E_w2, E_e2, ...
                                              fs_P2, fs_w2, fs_e2, ...
                                              A_eff_edges_xt2(time2-1:time2,:), ...
                                              b_dot_P2(time2-1:time2,:), ...
                                              b_dot_edges_SS2, ...
                                              b_dot_edges2(time2-1:time2,:), ...
                                              b_dot_P2, b_dot_edges2, ...
                                              Q_out_L_SS2, Q_out_R_SS2, ...
                                              S_at_GL2(time2-1:time2), ...
                                              Q_ext_P2(time2-1:time2,1), ...
                                              Q_ext_P2(time2-1:time2,2), ...  
                                              dt_P2(time2), t_P2, time2, ...
                                              deformation_only2, deformation_plus_sliding2, sliding_only2);
            
                                          
% fill matrices to save values
% ============================
     S_P2(time2,:)                = h_P2(time2,:) + B_P2; 
    
     h_dot2(time2,:)              = (h_P2(time2,:) - h_P2(time2-1,:)) / dt_P2(time2);
     
     S_dot2(time2,:)              = (S_P2(time2,:) - S_P2(time2-1,:)) / dt_P2(time2);
     
     [ dS_dx_w2, dS_dx_e2 ]       = get_gradient_values( S_P2(time2,:), x_P2, dx_P2 );
     dS_dx_edges_xt2(time2,:)     = [dS_dx_w2(1) dS_dx_e2];
     
     dS_dx_P_xt2(time2,:)         = qinterp1(x_edges2, dS_dx_edges_xt2(time2,:), x_P2)';
    
     
% FLUX
% ====
  flux_kin_P_xt2(time2,:) = calc_flux_kin( x_P2, x_P2, dx_P2, W_P2, ...
                                        h_dot2(time2,:), b_dot_P2(time2,:), Q_0_in2 );
                                                           
                                         
[ flux_w2, flux_e2 ] = get_edge_values_quadratic ( flux_kin_P_xt2(time2,:), ...
                                                  x_P2, x_w2, x_e2, ...
                                                  dx_P2, dx_w2, dx_e2 );
 
% % If spatial step is too large, may need to do this...
% if (flux_w(1) ~= Q_0_in)
%     flux_w(1)
%     flux_w(1) = Q_0_in;
%     disp(' Kinematic flux at downstream boundary does not equal Q_0_in, likely because spatial step too large')
%     disp(' Check at line 257 in calculate_surface')
%     disp(' ')
% end
% 
%   flux_edges_kin_xt(time,:) = [ flux_w(1) flux_e ];

 flux_edges_kin_xt2(time2,:) = [ Q_0_in2 flux_e2];   % Need to replace with Q_0_in, not extrapolated value.
                                                    

%  flux_edges_dyn_xt2(time2,:) = calc_flux_dyn( x_edges2, [h_w2(time2,1) h_e2(time2,:)], ...
%                                             dS_dx_edges_xt2(time2,:), ...
%                                             [E_w2(1) E_e2], [fs_w2(1) fs_e2], ...
%                                             [W_w2(1) W_e2], A_eff_edges_xt2(time2,:), ...
%                                             deformation_only2, deformation_plus_sliding2, sliding_only2);   
 
% flux_edges_dyn_xt2(1,1)    = Q_0_in2; 

 
 average_vel_estimate2 = abs(flux_edges_dyn_xt2(time2,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
                                        
                                        

rmpath( DIRECTORY_surf )    % remove path of MFILES_surface_main    



end   % if statement on time



z_grid2 = ( (z_hat_grid2' .* repmat(h_P2(time2,:), N_z2, 1)) + repmat(B_P2, N_z2, 1) )';


% ____________________________________________________   end loop on time     



% Save file with no parameter changes to compare against in figures:
if ((min_search_bed == 0) && (min_search_E == 0) && (min_search_fs == 0))
    
x_P_nochanges                  = x_P;
S_P_nochanges                  = S_P;
x_edges_nochanges              = x_edges;
average_vel_estimate_nochanges = average_vel_estimate;

save values_nochanges.mat x_P_nochanges S_P_nochanges x_edges_nochanges average_vel_estimate_nochanges

end




%    [ u, w, w_dzhat_dt, ...
%      u_bar_P, ...
%      u_bar_edges, phi, ...
%      psi, int_dphi_dx,  ...
%      A_eff, ...
%      wterm1, wterm2, ...
%      wterm3, wterm4 ] =  velocity_field ( x_grid', z_grid', z_hat_grid', T_use', ...
%                                           dS_dx_P_xt(time,:), ...
%                                           h_P(time,:), [h_w(time,1) h_e(time,:)], ...
%                                           dB_dx_P, b_dot_P(time,:), ...
%                                           S_dot(time,:), B_P, S_P(time,:), dz_hat, ...
%                                           A_0, x_P, x_w, x_e, dx_P, dx_w, dx_e, ...
%                                           flux_kin_P_xt(time,:), ...
%                                           flux_edges_kin_xt(time,:), ...
%                                           W_P, W_w, W_e );             
