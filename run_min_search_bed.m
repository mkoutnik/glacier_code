% run_min_search_bed


bed_vec         = [ 400:200:4000 ];
N_bed           = length( bed_vec );

length_run    = length(x_P);
save_bed        = NaN * ones(1, N_bed);
save_mismatch = NaN * ones(length_run, N_bed);

B_P_orig = B_P;

spliced_bed       = 151100;  % at this distance along flowband
spliced_bed_index = find(x_P <= spliced_bed, 1, 'last');  % Give one more grid point?


%for xpos = length_run:-1:spliced_bed_index   % loop over unknown bed positions
%for xpos = spliced_bed_index:length_run 
for xpos = length_run:-1:1 

disp(['Evaluating bed for x-position ',int2str(xpos), ' of total=',int2str(length_run)]);    
    
RMS_mismatch = NaN * ones( 1, N_bed );


%  Loop over B_P
 for i_bed = 1:N_bed

      B_P     = B_P_orig; % RESET
      bed_use = bed_vec(i_bed);

      B_P(xpos) = bed_use;
      
    % recalculate values on the edges.
      [B_w, B_e ] = get_edge_values_quadratic ...
                             ( B_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

   [ S_P(time,:), h_P(time,:), ...  
    dS_dx_P_xt(time,:),  ...
    dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
                                         W_P, W_w, W_e, b_dot_edges(time,:), ...
                                         E_w, E_e, fs_w, fs_e, ...
                                         Q_0_in, S_at_GL(time), ...
                                         A_eff_edges_xt(time,1:end-1), ...
                                         A_eff_edges_xt(time,2:end), ...
                                         flux_add_w, flux_add_e);
                                     
                         
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
  % residual = sqrt( mean( ( (abs(BDM_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
        
 % Only at xpos point?? 
   residual = sqrt( mean( ( (abs(BDM_velocity(xpos) - abs(surf_vel_estimate(xpos))/std_dev ).^2 ) ) ) );
   
   RMS_mismatch(i_bed) = residual;
   
                                   
 end  % loop over E values
 
 
   index     = find( RMS_mismatch == min(min(RMS_mismatch) ) );
   RMS_best  = RMS_mismatch(index);
   best_bed  = bed_vec(index)
   save_bed(xpos)  = best_bed;    
   save_mismatch(xpos, :) = RMS_mismatch;

% set best bed:
%  B_P(xpos) = best_bed;
 
end  % loop over xpositions

                 

% B_P_min = [B_P(1:spliced_bed_index-1) save_bed];
            
B_P_min = save_bed;

B_P = B_P_min;
% recalculate values on the edges.
      [B_w, B_e ] = get_edge_values_quadratic ...
                             ( B_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

figure(100)
plot(x_P, B_P_orig,'k')
hold on
plot(x_P, B_P_min, 'r--')
title('Bed')

% -------------------------------------------------------------------------