% script.           


global dt_impulse
global N_t_mesh
global all_IRF




% The impulse response function is the flux variation across the
% limited-domain boundaries in response to an impulsive perturbation in
% accumulation
% ----------------------------------------------------------------------



% cutoff response function when returned to steady state
% =======================================================
  test_stop_L = abs( h_P_xt(time_perturb:end,index_ldm_L) - ...
                     h_P_xt(time_perturb-1:end-1,index_ldm_L) );  
                 
  test_stop_R = abs( h_P_xt(time_perturb:end,index_ldm_R) - ...
                     h_P_xt(time_perturb-1:end-1,index_ldm_R) );  
  
  start_use   = 50;   % problems with cutoff too early if use time_perturb as
                      % the start value for impulsive perturbation (in
                      % generating the divide IRF). try to avoid confusion.
                 
  time_stop_L = find(test_stop_L(start_use:end) <= ss_return_value, 1 ) + start_use;
  time_stop_R = find(test_stop_R(start_use:end) <= ss_return_value, 1 ) + start_use;
 
  
  
if ( (isempty(time_stop_L)) || (isempty(time_stop_R)) )
    disp(' ')
    disp('USING ALL CALCULATED VALUES OF h1 in RESPONSE_FUNCTION_CALC.m')
    stop;
end




 
% flux response from numerical value
% (dynamic flux calculation)
% ==========================
  h1_ldm_L            = h_P_xt(time_perturb:time_stop_L, index_ldm_L);
  h1_ldm_R            = h_P_xt(time_perturb:time_stop_R, index_ldm_R);
  
  q1_ldm_L            = -flux_edges_dyn_xt(time_perturb+1:time_stop_L, index_ldm_L);
  
  q1_ldm_normalized_L = (q1_ldm_L - min(q1_ldm_L)) / ...
                         max(q1_ldm_L - min(q1_ldm_L)); % integral from 0-tau = 1
  
  q1_ldm_R            = flux_edges_dyn_xt(time_perturb+1:time_stop_R, index_ldm_R);
  
  q1_ldm_normalized_R = (q1_ldm_R - min(q1_ldm_R)) / ...
                         max(q1_ldm_R - min(q1_ldm_R)); % integral from 0-tau = 1 
                     
   
 
                     
% time until returns to steady state
% ==================================
  tau_L = (time_stop_L-time_perturb) * dt_impulse;
  tau_R = (time_stop_R-time_perturb) * dt_impulse;
  
  

% timestep for response calculation
% =================================
  t_grid_imp_L = time_vec(time_perturb-1:time_stop_L);
  t_grid_imp_R = time_vec(time_perturb-1:time_stop_R);

  
% compare to volume response time, H_max/bdot_term
% ================================================
  tau_volume_L = max(h_P) / abs(ablation_rate_L);
  tau_volume_R = max(h_P) / abs(ablation_rate_R);
  
 

    

% dt for q1_ldm, and dt for transient run must be the same.
  % ---------------------
  if (dt_impulse ~= dt_P)   % need to interpolate onto mesh
  % ---------------------

% for the left:
% =============
   t_P_interp_L_temp = [ 0: dt_P: t_grid_imp_L(end-2)]; % full grid for response time 
                
%      if (length(t_P_interp_L) > N_t_mesh)        % only keep relevant values
%          t_P_interp_L = [ 0:dt_P:t_P_interp_L(N_t_mesh) ]; 
%      end
   
     q1_on_dt_use_L = interp1(t_grid_imp_L-t_grid_imp_L(1), q1_ldm_L, ...
                              t_P_interp_L_temp);
     % there is no extrapolation here, so if length(t_P_interp_L) > length(t_grid_imp_L),
     % the value will be NaN... change to zero later!                     
                          
     h1_on_dt_use_L = interp1(t_grid_imp_L-t_grid_imp_L(1), h1_ldm_L, ...
                              t_P_interp_L_temp);
     
     
% for the right:
% ==============
   t_P_interp_R_temp = [ 0: dt_P: t_grid_imp_R(end-2)]; % full grid for response time
                
%      if (length(t_P_interp_R) > N_t_mesh)          % only keep relevant values
%          t_P_interp_R = [ 0:dt_P:t_P_interp_R(N_t_mesh) ];
%      end
  
     q1_on_dt_use_R = interp1(t_grid_imp_R-t_grid_imp_R(1), q1_ldm_R, ...
                              t_P_interp_R_temp);    
     % there is no extrapolation here, so if length(t_P_interp_R) > length(t_grid_imp_R),
     % the value will be NaN... change to zero later!
     
     h1_on_dt_use_R = interp1(t_grid_imp_R-t_grid_imp_R(1), h1_ldm_R, ...
                              t_P_interp_R_temp);    
     
% -------                       
  else
% -------                       
    
   q1_on_dt_use_L = q1_ldm_L';
   q1_on_dt_use_R = q1_ldm_R';
   
   h1_on_dt_use_L = h1_ldm_L';
   h1_on_dt_use_R = h1_ldm_R';
   

   t_P_interp_L_temp = t_grid_imp_L;
   t_P_interp_R_temp = t_grid_imp_R;
   
   
  end     % end on dt_imp check.   
  
  
  
% But...
% ======
% Replace NaN. 
  if (sum(isnan(q1_on_dt_use_L)) > 0)
      disp('NaN in response_function_calc.m')
      stop;
%        q1_on_dt_use_L(isnan(q1_on_dt_use_L) == 1) = 0;
  end
  
  if (sum(isnan(q1_on_dt_use_R)) > 0)
      disp('NaN in response_function_calc.m')
      stop;
      % q1_on_dt_use_R(isnan(q1_on_dt_use_R) == 1) = 0;
  end
               
  

  
% the actual response function
% =============================
% q1_on_dt_use_L is negative, but IRF should be positive because volume change
% carries the sign. IRF just scales the volume change over time.

  q1_on_dt_L_temp = (abs(q1_on_dt_use_L) - min(abs(q1_on_dt_use_L))) / ...
                    sum( (abs(q1_on_dt_use_L)) - min(abs(q1_on_dt_use_L)) );  % scale so that integral = 1
  
  q1_on_dt_R_temp = (q1_on_dt_use_R - min(q1_on_dt_use_R)) / ...
                    sum((q1_on_dt_use_R - min(q1_on_dt_use_R)));  % scale so that integral = 1
    
                
 % this is only used for plotting, and comparison.               
  h1_on_dt_L_temp = (abs(h1_on_dt_use_L) - min(abs(h1_on_dt_use_L))) / ...
                    sum( (abs(h1_on_dt_use_L)) - min(abs(h1_on_dt_use_L)) );  % scale so that integral = 1
  
  h1_on_dt_R_temp = (h1_on_dt_use_R - min(h1_on_dt_use_R)) / ...
                    sum((h1_on_dt_use_R - min(h1_on_dt_use_R)));  % scale so that integral = 1
    
                



% size of q1 must be size N_t_mesh, so fill with zeros:
if (length(q1_on_dt_L_temp) < N_t_mesh)
    q1_on_dt_L   = [ q1_on_dt_L_temp zeros(1, N_t_mesh-length(q1_on_dt_L_temp)) ]';
    h1_on_dt_L   = [ h1_on_dt_L_temp zeros(1, N_t_mesh-length(h1_on_dt_L_temp)) ]';
    t_P_interp_L = [ t_P_interp_L_temp NaN(1, N_t_mesh-length(t_P_interp_L_temp)) ]';
elseif (length(q1_on_dt_L_temp) > N_t_mesh)
    q1_on_dt_L   = q1_on_dt_L_temp(1:N_t_mesh);
    h1_on_dt_L   = h1_on_dt_L_temp(1:N_t_mesh);
    t_P_interp_L = t_P_interp_L_temp(1:N_t_mesh);
end


if (length(q1_on_dt_R_temp) < N_t_mesh)
    q1_on_dt_R   = [ q1_on_dt_R_temp zeros(1, N_t_mesh-length(q1_on_dt_R_temp)) ]';
    h1_on_dt_R   = [ h1_on_dt_R_temp zeros(1, N_t_mesh-length(h1_on_dt_R_temp)) ]';
    t_P_interp_R = [ t_P_interp_R_temp NaN(1, N_t_mesh-length(t_P_interp_R_temp)) ]';
elseif (length(q1_on_dt_R_temp) > N_t_mesh)
    q1_on_dt_R   = q1_on_dt_R_temp(1:N_t_mesh);
    h1_on_dt_R   = h1_on_dt_R_temp(1:N_t_mesh);
    t_P_interp_R = t_P_interp_R_temp(1:N_t_mesh);
end
           








if (all_IRF == 1)      % to scale ice-divide migration perturbations
    
 save all_impulse_response_function_values_div.mat ...
        index_ldm_L index_ldm_R test_stop_L test_stop_R ...
        time_stop_L time_stop_R ...
        h1_ldm_L h1_ldm_R ...
        q1_ldm_L q1_ldm_R q1_ldm_normalized_L ...
        q1_ldm_normalized_R tau_L tau_R t_grid_imp_L t_grid_imp_R ...
        tau_volume_L tau_volume_R    
    
    
    
 % this gets loaded later in convolution.m           
  save q1_on_dt_div.mat q1_on_dt_L h1_on_dt_L h1_on_dt_R tau_L q1_on_dt_R ...
                        tau_R t_P_interp_R t_P_interp_L    
                      
                    
  
elseif (all_IRF == 2)   % to scale accumulation perturbations   
   
  save all_impulse_response_function_values_vol.mat ...
       index_ldm_L index_ldm_R test_stop_L test_stop_R ...
        time_stop_L time_stop_R ...
        h1_ldm_L h1_ldm_R ...
        q1_ldm_L q1_ldm_R q1_ldm_normalized_L ...
        q1_ldm_normalized_R tau_L tau_R t_grid_imp_L t_grid_imp_R ...
        tau_volume_L tau_volume_R   
    
  % this gets loaded later in convolution.m           
  save q1_on_dt_vol.mat q1_on_dt_L h1_on_dt_L h1_on_dt_R ...
                        tau_L q1_on_dt_R tau_R t_P_interp_R t_P_interp_L
  
end


  
% check to make sure IRF is reasonable:
  if ( (abs(q1_on_dt_L(2)-q1_on_dt_L(1)) < 1e-6) || (abs(q1_on_dt_R(2)-q1_on_dt_R(1)) < 1e-6) )
      disp('Response function incorrect. CHECK response_function_calc.m')
      stop;
  end
  
  
  
  
  
  
  
  
  
  