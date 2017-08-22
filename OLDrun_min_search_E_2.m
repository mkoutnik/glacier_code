% run_min_search_E

% See if E patterned after the observed surface velocity is an
% improvement... not really... compromise surface velocity fit for surface
% topography



E_vec         = [ 0.1 1 2 3 5 6 7 8 9 10 20 50 ];
N_E           = length( E_vec );

length_run    = length(x_P);

% E_P_orig = E_P;



BDM_velocity        = interp1(measures_centerline_distance, measures_flowspeed, x_P);
E_P_orig            = BDM_velocity/max(BDM_velocity);


% % E_test_values = NaN * ones(N_E*length_run, length_run);
% % 
% % % setup field of possible E values
% % for ii = 1:N_E
% %     for jj = 1:length_run
% %  
% %         E_test_values(ii, jj) = E_vec(ii);
% %     
% %     end
% % end


% for xpos = length_run:-1:1   % loop over all x positions -- downstream
% % for xpos = 1:length_run   % loop over all x positions -- upstream
    
% disp(['Evaluating E for x-position ',int2str(xpos), ' of total=',int2str(length_run)]);

RMS_mismatch = NaN * ones( 1, N_E );


%  Loop over E_x -- this time it is a factor to E(x) profile from surf vel
 for i_E = 1:N_E

      E_use = E_vec(i_E);

     % E_P(xpos) = E_use;
       E_P = E_P_orig * E_use;  % multiply function by E factor 
     
    % recalculate values on the edges.
      [ E_w, E_e ] = get_edge_values_quadratic ...
                               ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

                           
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

 % Surface profile:
  % residual = sqrt( mean( ( (abs(S_modern) - abs(S_P(1,:)))/std_dev ).^2 ) );           

 % Surface velocity:
  % residual = sqrt( mean( ( (abs(BDM_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
  % %  residual = abs(BDM_velocity(xpos) - surf_vel_estimate(xpos));
  
  
  
  residual = sqrt( mean( ( (abs(BDM_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) ) + ...
             sqrt( mean( ( (abs(S_modern) - abs(S_P(1,:)))/std_dev ).^2 ) );  
  
  
   RMS_mismatch(i_E) = residual;
   
                                   
 end  % loop over E values
 
 
   index     = find( RMS_mismatch == min(min(RMS_mismatch) ) );
   
   if (length(index)>1)  % Multiple values of E give same mismatch
       disp('Multiple values of E give same mismatch value -- set E=1')
       index = find(E_vec == 1);
   end
   
   RMS_best  = RMS_mismatch(index)
   best_E    = E_vec(index)

   
 % Keep previous values:  
%   E_P(xpos) = best_E;
  
 % % OR, reset E:
 %  E_P   = E_P_orig;
     
 
% end  % loop over xpositions



       E_P = best_E * E_P_orig;
       
     % recalculate values on the edges.
       [ E_w, E_e ] = get_edge_values_quadratic ...
                                ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );



% % -------------------------------------------------------------------------
% 
% % run_min_search_E
% 
% E_vec         = [ 0.1 1 2 3 4 5 6 7 8 9 10 ];
% N_E           = length( E_vec );
% 
% length_run    = length(x_P);
% save_E        = [];
% save_mismatch = NaN * ones(length_run, N_E);
% 
% E_P_orig = E_P;
% 
% 
%   for xpos = length_run:-1:1   % loop over all x positions
% % for xpos = 1:length_run   % loop over all x positions
%     
% disp(['Evaluating for x-position ',int2str(xpos), ' out of total=',int2str(length_run)]);
% 
% RMS_mismatch = NaN * ones( 1, N_E );
% 
% 
% %  Loop over E_x
%  for i_E = 1:N_E
% 
%       E_use = E_vec(i_E);
% 
%       E_P(xpos) = E_use;
%       
%     % recalculate values on the edges.
%       [ E_w, E_e ] = get_edge_values_quadratic ...
%                                ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
% 
%                            
%   [ S_P(time,:), h_P(time,:), ...  
%     dS_dx_P_xt(time,:),  ...
%     dS_dx_edges_xt(time,:) ] = calc_h_0( x_P, x_w, x_e, dx_P, B_P, B_w, B_e, ...
%                                          W_P, W_w, W_e, b_dot_edges(time,:), ...
%                                          E_w, E_e, fs_w, fs_e, ...
%                                          Q_0_in, S_at_GL(time), ...
%                                          A_eff_edges_xt(time,1:end-1), ...
%                                          A_eff_edges_xt(time,2:end) );
%                                      
%                          
% % edge values of ice thickness
% % ============================
%   [ h_w(time,:), ...
%     h_e(time,:) ] = get_edge_values_quadratic( h_P(time,:), x_P, x_w, x_e, ...
%                                                dx_P, dx_w, dx_e );
%               
% % Find the flux
% % ==============
%  flux_kin_P_xt(time,:) = calc_flux_kin( x_P, x_P, dx_P, W_P, ...
%                                         h_dot(time,:), b_dot_P(time,:), Q_0_in );
%                                                                                              
% [ flux_w, flux_e ] = get_edge_values_quadratic ( flux_kin_P_xt(time,:), ...
%                                                  x_P, x_w, x_e, ...
%                                                  dx_P, dx_w, dx_e );
%    
%  flux_edges_kin_xt(time,:) = [ Q_0_in flux_e];   % Need to replace with Q_0_in, not extrapolated value.
%                                                     
%  flux_edges_dyn_xt(time,:) = calc_flux_dyn( x_edges, [h_w(time,1) h_e(time,:)], ...
%                                             dS_dx_edges_xt(time,:), ...
%                                             [E_w(1) E_e], [fs_w(1) fs_e], ...
%                                             [W_w(1) W_e], A_eff_edges_xt(time,:) );   
%  flux_edges_dyn_xt(1,1)    = Q_0_in; 
%                                         
%            
%  
%  % Velocity
%  % ========
%  surf_vel_estimate = (5/4) * abs(flux_edges_dyn_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
%  
%  BDM_velocity      = interp1(measures_centerline_distance, measures_flowspeed, [x_w(1) x_e]);
%  
%  
%  % Residual
%  % ========
%     std_dev = 1; 
% 
%  % Surface profile:
%   residual = sqrt( mean( ( (abs(S_modern) - abs(S_P(1,:)))/std_dev ).^2 ) );           
% 
%  % Surface velocity:
%  %  residual = sqrt( mean( ( (abs(BDM_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
%         
%    RMS_mismatch(i_E) = residual;
%    
%                                    
%  end  % loop over E values
%  
%  
%    index     = find( RMS_mismatch == min(min(RMS_mismatch) ) );
%    RMS_best  = RMS_mismatch(index);
%    best_E    = E_vec(index);
%    save_E    = [ save_E best_E ];    
%    save_mismatch(xpos, :) = RMS_mismatch;
% 
% % reset E:
%   E_P(xpos) = best_E;
%  
% end  % loop over xpositions
% 
% 
% 
% 
% 
% % -------------------------------------------------------------------------
