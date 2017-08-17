% run_min_search_E

E_vec         = [ 0.01 0.1 0.5 1 2 4 5 6 8 10 20 50 100 ];
N_E           = length( E_vec );


% CHANGE iterations = 1 or = 2 here!!
% Might be working to run through both together but don't have file saving
% setup and don't have it double checked. 


for iterations = 1   % Minimize for both Darwin (1) and Hatherton (2)
% for iterations = 2    
    
    if (iterations == 1)
        disp('Running for Darwin...')
        length_run      = length(x_P);
        save_mismatch   = NaN * ones(length_run, N_E);
        save_E          = NaN * ones(1, N_E);
        E_P_orig        = E_P;
        Darwin_velocity = interp1(Darwin_measures_centerline_distance, Darwin_measures_flowspeed, [x_w(1) x_e]);
        Darwin_surface  = S_modern;

        
    

 for xpos = 1:length_run   % loop over all x positions -- upstream
    
disp(['Evaluating E for x-position ',int2str(xpos), ' of total=',int2str(length_run)]);

RMS_mismatch = NaN * ones( 1, N_E );


%  Loop over E_x
 for i_E = 1:N_E
     
      E_use = E_vec(i_E);
      E_P(xpos) = E_use;
      
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
                                         flux_add_w, flux_add_e, ...
                                         deformation_only, deformation_plus_sliding, sliding_only);
                                                         
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
                                            [W_w(1) W_e], A_eff_edges_xt(time,:), ...
                                            deformation_only, deformation_plus_sliding, sliding_only);   
 flux_edges_dyn_xt(1,1)    = Q_0_in; 
           
 
 % USING KIN FLUX HERE!
           
 % Velocity
 % ========
 surf_vel_estimate = (5/4) * abs(flux_edges_kin_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
 average_vel_estimate = (5/4) * abs(flux_edges_kin_xt(1,:)) ./ ([W_w(1) W_e] .* [h_w(1,1) h_e(1,:)]);
 
 
 % Residual
 % ========
    std_dev = 1; 

%  % Surface profile:
%   residual = sqrt( mean( ( (abs(Darwin_surface) - abs(S_P(1,:)))/std_dev ).^2 ) );           
%    residual = sqrt( mean( ( (abs(Darwin_surface(xpos)) - abs(S_P(1,xpos)))/std_dev ).^2 ) );  
         

 % % Surface velocity:
  % residual = sqrt( mean( ( (abs(Darwin_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) );  
  %  residual = abs(Darwin_velocity(xpos) - surf_vel_estimate(xpos));  % local mismatch?
 
 % Both:
  %  residual = sqrt( mean( ( (abs(Darwin_velocity) - abs(surf_vel_estimate))/std_dev ).^2 ) ) + ...
  %             sqrt( mean( ( (abs(Darwin_surface) - abs(S_P(1,:)))/std_dev ).^2 ) );    
    residual = sqrt( mean( ( (abs(Darwin_velocity(xpos)) - abs(surf_vel_estimate(xpos)))/std_dev ).^2 ) ) + ...
               sqrt( mean( ( (abs(Darwin_surface(xpos)) - abs(S_P(1,xpos)))/std_dev ).^2 ) );   
      
  
if (isreal(S_P(1,:)) == 0)
    residual = NaN;   % Don't pick from this value if gives imaginary.
end
  
   RMS_mismatch(i_E) = residual;
   
                                   
 end  % loop over E values
 
 
 
   index     = find( RMS_mismatch == min(min(RMS_mismatch) ) );
   
   if (length(index)>1)  % Multiple values of E give same mismatch
       disp('Multiple values of E give same mismatch value -- set E=1')
       index = find(E_vec == 1);
   end
   
   RMS_best  = RMS_mismatch(index);
   best_E    = E_vec(index);
   save_E(xpos) = best_E;    
   save_mismatch(xpos, :) = RMS_mismatch;

   
 % Keep previous values:  
   E_P(xpos) = best_E;
  
 % % OR, reset E:
 %  E_P   = E_P_orig;
     
 
end  % loop over xpositions


      E_P = save_E;
      
    % recalculate values on the edges.
      [ E_w, E_e ] = get_edge_values_quadratic ...
                               ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

                        
%       calculate_surface_Darwin;
       
% % RESET Hatherton from Darwin calculation
% S_in_global_Hatherton = h_w(1, index_xpos_Hatherton_edges) + B_w(1, index_xpos_Hatherton_edges);
% S_0_in2               = S_in_global_Hatherton;
% S_at_GL2              = ones(size(t_P2)) * S_0_in2;  
%       
%      
% 
%     calculate_surface_Hatherton;


   
   


 % ------------------------------------------------------------------------                          
  elseif (iterations == 2)
        disp('Running for Hatherton...')
        length_run2      = length(x_P2);
        save_mismatch2   = NaN * ones(length_run2, N_E);
        save_E2          = NaN * ones(1, N_E);
        E_P_orig2        = E_P2;
        Hat_velocity    = interp1( Hat_measures_centerline_distance, Hat_measures_flowspeed, x_edges2, 'linear', 'extrap');
        Hat_surface     = S_modern2;

                            

 for xpos2 = 1:length_run2   % loop over all x positions -- upstream
    
disp(['Evaluating E for x-position ',int2str(xpos2), ' of total=',int2str(length_run2)]);

RMS_mismatch2 = NaN * ones( 1, N_E );


%  Loop over E_x
 for i_E = 1:N_E
     
      E_use = E_vec(i_E);
      E_P2(xpos2) = E_use;
      
    % recalculate values on the edges.
      [ E_w2, E_e2 ] = get_edge_values_quadratic ...
                               ( E_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );
                         
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
    
  
 if (isreal(S_P2(time2,:)) == 0)
   disp('Imaginary in surface values')
 end

                                     
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
                                                    
 flux_edges_dyn_xt2(time2,:) = calc_flux_dyn( x_edges2, [h_w2(time2,1) h_e2(time2,:)], ...
                                            dS_dx_edges_xt2(time2,:), ...
                                            [E_w2(1) E_e2], [fs_w2(1) fs_e2], ...
                                            [W_w2(1) W_e2], A_eff_edges_xt2(time2,:), ...
                                            deformation_only2, deformation_plus_sliding2, sliding_only2);   
                                    
                                     
 flux_edges_dyn_xt2(1,1)    = Q_0_in; 
                                        
         
 % USING KIN FLUX HERE!
 
 % Velocity
 % ========
 surf_vel_estimate2 = (5/4) * abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
 average_vel_estimate2 = abs(flux_edges_kin_xt2(1,:)) ./ ([W_w2(1) W_e2] .* [h_w2(1,1) h_e2(1,:)]);
 
 
 % Residual
 % ========
    std_dev = 1; 

 % % Surface profile:
 %  residual2 = sqrt( mean( ( (abs(Hat_surface) - abs(S_P2(1,:)))/std_dev ).^2 ) );           

 % Surface velocity:
  % residual2 = sqrt( mean( ( (abs(Hat_velocity) - abs(surf_vel_estimate2))/std_dev ).^2 ) );  
  %  residual2 = abs(Hat_velocity(xpos2) - surf_vel_estimate2(xpos2));  % local mismatch?
 
 % % Both:
 %   residual2 = sqrt( mean( ( (abs(Hat_velocity) - abs(surf_vel_estimate2))/std_dev ).^2 ) ) + ...
 %            sqrt( mean( ( (abs(Hat_surface) - abs(S_P2(1,:)))/std_dev ).^2 ) );    
   residual2 = sqrt( mean( ( (abs(Hat_velocity(xpos2)) - abs(surf_vel_estimate2(xpos2)))/std_dev ).^2 ) ) + ...
              sqrt( mean( ( (abs(Hat_surface(xpos2)) - abs(S_P2(1,xpos2)))/std_dev ).^2 ) );   
          
   
  
   RMS_mismatch2(i_E) = residual2;
   
                                   
 end  % loop over E values
 
 
 
   index2     = find( RMS_mismatch2 == min(min(RMS_mismatch2) ) );
   
   if (length(index2)>1)  % Multiple values of E give same mismatch
       disp('Multiple values of E give same mismatch value -- set E=1')
       index2 = find(E_vec == 1);
   end
   
   RMS_best2  = RMS_mismatch2(index2);
   best_E2    = E_vec(index2);
   save_E2(xpos2) = best_E2;    
   save_mismatch2(xpos2, :) = RMS_mismatch2;

   
 % Keep previous values:  
   E_P2(xpos2) = best_E2;
  
 % % OR, reset E:
 %  E_P   = E_P_orig;
     
 
end  % loop over xpositions


      E_P2 = save_E2;
      
    % recalculate values on the edges.
      [ E_w2, E_e2 ] = get_edge_values_quadratic ...
                               ( E_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );            
        
        
        
        
% % Then, setup flux added at 17 km from Hatherton...
% flux_add_P(index_xpos_Hatherton_P) = flux_edges_dyn_xt2(time, 1) / W_P2(1);
% [ flux_add_w, ...
%   flux_add_e ] = get_edge_values_quadratic ...
%                                   ( flux_add_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
% 
%                               
%    calculate_surface_Darwin;
%       disp('Recalculating for Darwin')

                           
    end                        
 
end
 
 
 