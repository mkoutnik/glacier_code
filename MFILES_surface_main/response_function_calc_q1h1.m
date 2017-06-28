global n
global dt_impulse
global N_t_mesh
global all_IRF



% ------------------------------------
% calculate impulse response function
% ------------------------------------

% C_0_ldm -- kinematic wave velocity
% D_0_ldm -- wave diffusion coefficient
% definitions from e.g. Hooke (2005) pg. 371, 375



% values for the limited domain
% ==============================
  [ h_w, h_e ]         = get_edge_values_quadratic( h_P, x_P, x_w, x_e, ...
                                                    dx_P, dx_w, dx_e );
                                        
  [ B_w, B_e ]         = get_edge_values_quadratic( B_P, x_P, x_w, x_e, ...
                                                    dx_P, dx_w, dx_e );
                                                
                                                
  [ dh_dx_w, dh_dx_e ] = get_gradient_values( h_P, x_P, dx_P );

  
  [ dB_dx_w, dB_dx_e ] = get_gradient_values( B_P, x_P, dx_P );
 
  
% alpha = dB_dx - dh_dx   (Hooke pg. 368)  
  alpha_w = dB_dx_w - dh_dx_w;
  alpha_e = dB_dx_e - dh_dx_e;
  
  
% kinematic wave velocity:  
  C_0_all = (n+2) * u_bar_edges;
  
% kinematic wave diffusion coefficient:  
  D_0_all = (n * u_bar_edges .* [h_w(1) h_e]) ./ [alpha_w(1) alpha_e];
  
  
% these are on the limited domain only:  
  C_0_ldm_L = (n+2) * u_bar_edges(1);   
  D_0_ldm_L = (n * u_bar_edges(1) * h_w(1)) / ...
               alpha_w(1);
  
  C_0_ldm_R = (n+2) * u_bar_edges(end);   
  D_0_ldm_R = (n * u_bar_edges(end) * h_e(end)) / ...
               alpha_e(end);
  
% check the signs:
C_0_ldm_L = -abs(C_0_ldm_L);
D_0_ldm_L = -abs(D_0_ldm_L);
C_0_ldm_R = abs(C_0_ldm_R);
D_0_ldm_R = -abs(D_0_ldm_R);



% % Or, use from updated values?
% % ----------------------------
%                                         
% [ dh_dx_xt_w, dh_dx_xt_e ] = get_gradient_values( h_P_xt(1,:), full_x_P, full_dx_P );
% 
% C_0_ldm_L = (n+2) * (flux_edges_xt(1,index_ldm_L)/h_e_xt(1,index_ldm_L));
% D_0_ldm_L = (-n * (flux_edges_xt(1,index_ldm_L))/ (dh_dx_xt_e(index_ldm_L)) );
%
% C_0_ldm_R = (n+2) * (flux_edges_xt(1,index_ldm_R)/h_e_xt(1,index_ldm_R));
% D_0_ldm_R = (-n * (flux_edges_xt(1,index_ldm_R))/ (dh_dx_xt_e(index_ldm_R)) );

                            

% cutoff response function when returned to steady state
% =======================================================
  test_stop_L = abs( h_P_xt(time_perturb+1:end,index_ldm_L) - ...
                     h_P_xt(time_perturb:end-1,index_ldm_L) );  
                 
  test_stop_R = abs( h_P_xt(time_perturb+1:end,index_ldm_R) - ...
                     h_P_xt(time_perturb:end-1,index_ldm_R) );  
  
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
 
 
 
  
% h1(index_ldm, t) -- changes across ldm boundary
% ===============================================


% for the left:
% -------------
  h1_ldm_all_L = h_P_xt(:,index_ldm_L); 
  h1_ldm_L     = (h1_ldm_all_L(time_perturb:time_stop_L)) - h1_ldm_all_L(1); 
  
 
  
% % recalculate time_stop to be equal 2*tau (e-folding time)                     
%   h1_response_L = (h1_ldm_L - min(h1_ldm_L)) / ...
%                  sum((h1_ldm_L - min(h1_ldm_L)));  % scale so that integral = 1
%                                 
%   efolding_index_L = (max(find(cumsum(h1_response_L) <= (1-(1/2.71828)))));
%   time_stop_L = 2 * efolding_index_L;
%   
%   h1_ldm_L     = h1_ldm_all_L(time_perturb:time_stop_L); 
  
  h1_ldm_normalized_L = (h1_ldm_L - min(h1_ldm_L)) / ...
                         max(h1_ldm_L - min(h1_ldm_L)); % integral from 0-tau = 1
                     
                     
% for the right:
% --------------                     
  h1_ldm_all_R = h_P_xt(:,index_ldm_R);
  h1_ldm_R     = h1_ldm_all_R(time_perturb:time_stop_R) - h1_ldm_all_R(1); 
  
  
%   % recalculate time_stop to be equal 2*tau (e-folding time)                     
%   h1_response_R = (h1_ldm_R - min(h1_ldm_R)) / ...
%                  sum((h1_ldm_R - min(h1_ldm_R)));  % scale so that integral = 1
%                                 
%   efolding_index_R = (max(find(cumsum(h1_response_R) <= (1-(1/2.71828)))));
%   time_stop_R = 2 * efolding_index_R;
%   
%   h1_ldm_R     = h1_ldm_all_R(time_perturb:time_stop_R); 
  
  h1_ldm_normalized_R = (h1_ldm_R - min(h1_ldm_R)) / ...
                         max(h1_ldm_R - min(h1_ldm_R)); % integral from 0-tau = 1
  
                     

% -----------------------------------
% flux response at LDM boundary, q1
% -----------------------------------

% slope, alpha (dh/dx)
% ====================
  h_diff_L     = (h_P_xt(time_perturb:time_stop_L,index_ldm_L) - h_P_xt(1,index_ldm_L)) - ...
                 (h_P_xt(time_perturb:time_stop_L,index_ldm_L+1) - h_P_xt(1,index_ldm_L+1));
  x_diff_L     = full_x_P(index_ldm_L+1) - full_x_P(index_ldm_L);
  
  
  alpha1_ldm_L = abs( h_diff_L ./ x_diff_L );   % positive on the left
  

  h_diff_R     = (h_P_xt(time_perturb:time_stop_R,index_ldm_R) - h_P_xt(1,index_ldm_R)) - ...
                 (h_P_xt(time_perturb:time_stop_R,index_ldm_R+1) - h_P_xt(1,index_ldm_R+1));
  x_diff_R     = full_x_P(index_ldm_R+1) - full_x_P(index_ldm_R);
  
  alpha1_ldm_R = -abs( h_diff_R ./ x_diff_R );   % negative on the right
  

  
% only keep the number of values that need to be calculated.
% flux response function, q1 at LDM boundary
% ==========================================


  q1_ldm_L            = abs( (C_0_ldm_L * h1_ldm_L) + (D_0_ldm_L * alpha1_ldm_L) );
  q1_ldm_normalized_L = (q1_ldm_L - min(q1_ldm_L)) / ...
                         max(q1_ldm_L - min(q1_ldm_L)); % integral from 0-tau = 1
  
  q1_ldm_R            = abs(( C_0_ldm_R * h1_ldm_R) + (D_0_ldm_R * alpha1_ldm_R) );
  q1_ldm_normalized_R = (q1_ldm_R - min(q1_ldm_R)) / ...
                         max(q1_ldm_R - min(q1_ldm_R)); % integral from 0-tau = 1 
 

 % on the left side of the divide:
 %   ubar is negative; alpha is positive (c0 and D0 are negative)
 
 % on the right side of the divide:
 %   ubar is positive; alpha is negative  (c0 is positive; D0 is negative)
  
                     
 
                     
% time until returns to steady state
% ==================================
  tau_L = (time_stop_L-time_perturb) * dt_impulse;
  tau_R = (time_stop_R-time_perturb) * dt_impulse;
  
  

% timestep for response calculation
% =================================
  t_grid_imp_L = time_vec(time_perturb:time_stop_L);
  t_grid_imp_R = time_vec(time_perturb:time_stop_R);

  
% compare to volume response time, H_max/bdot_term
% ================================================
  tau_volume_L = max(h_w) / abs(ablation_rate_L);
  tau_volume_R = max(h_w) / abs(ablation_rate_R);
  
 

    
   
% Setup response function that will be used
% =========================================

% dt for q1_ldm, and dt for transient run must be the same
  % ---------------------
  if (dt_impulse ~= dt_P)   % need to interpolate onto mesh
  % ---------------------

% first for the left:
% ===================
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
     
     
% then for the right:
% ===================
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
  
  q1_on_dt_R_temp = ((q1_on_dt_use_R) - min((q1_on_dt_use_R))) / ...
                    sum( ((q1_on_dt_use_R)) - min((q1_on_dt_use_R)) );  % scale so that integral = 1
    
                
 % this is only used for plotting, and comparison.               
  h1_on_dt_L_temp = (abs(h1_on_dt_use_L) - min(abs(h1_on_dt_use_L))) / ...
                    sum( (abs(h1_on_dt_use_L)) - min(abs(h1_on_dt_use_L)) );  % scale so that integral = 1
  
  h1_on_dt_R_temp = (h1_on_dt_use_R - min(h1_on_dt_use_R)) / ...
                    sum((h1_on_dt_use_R - min(h1_on_dt_use_R)));  % scale so that integral = 1
    
                
           

% size of q1 must be N_t_mesh, so fill with zeros:
if (length(q1_on_dt_L_temp) < N_t_mesh)
    q1_on_dt_L   = [ q1_on_dt_L_temp zeros(1, N_t_mesh-length(q1_on_dt_L_temp)) ]';
    h1_on_dt_L   = [ h1_on_dt_L_temp zeros(1, N_t_mesh-length(h1_on_dt_L_temp)) ]';
    t_P_interp_L = [ t_P_interp_L_temp NaN(1, N_t_mesh-length(t_P_interp_L_temp)) ]';
elseif (length(q1_on_dt_L_temp) >= N_t_mesh)
    q1_on_dt_L   = q1_on_dt_L_temp(1:N_t_mesh);
    h1_on_dt_L   = h1_on_dt_L_temp(1:N_t_mesh);
    t_P_interp_L = t_P_interp_L_temp(1:N_t_mesh);
end


if (length(q1_on_dt_R_temp) < N_t_mesh)
    q1_on_dt_R   = [ q1_on_dt_R_temp zeros(1, N_t_mesh-length(q1_on_dt_R_temp)) ]';
    h1_on_dt_R   = [ h1_on_dt_R_temp zeros(1, N_t_mesh-length(h1_on_dt_R_temp)) ]';
    t_P_interp_R = [ t_P_interp_R_temp NaN(1, N_t_mesh-length(t_P_interp_R_temp)) ]';
elseif (length(q1_on_dt_R_temp) >= N_t_mesh)
    q1_on_dt_R   = q1_on_dt_R_temp(1:N_t_mesh);
    h1_on_dt_R   = h1_on_dt_R_temp(1:N_t_mesh);
    t_P_interp_R = t_P_interp_R_temp(1:N_t_mesh);
end
           



if (all_IRF == 1)      % to scale ice-divide migration perturbations
    
 save all_impulse_response_function_values_div.mat ...
        C_0_all D_0_all index_ldm_L index_ldm_R C_0_ldm_L D_0_ldm_L ...
        C_0_ldm_R D_0_ldm_R test_stop_L test_stop_R time_stop_L time_stop_R ...
        h1_ldm_all_L h1_ldm_L h1_ldm_normalized_L h1_ldm_all_R h1_ldm_R ...
        h1_ldm_normalized_R q1_ldm_L q1_ldm_R q1_ldm_normalized_L ...
        q1_ldm_normalized_R tau_L tau_R t_grid_imp_L t_grid_imp_R ...
        tau_volume_L tau_volume_R alpha1_ldm_L alpha1_ldm_R   
    
    
    
 % this gets loaded later in convolution.m           
  save q1_on_dt_div.mat q1_on_dt_L h1_on_dt_L h1_on_dt_R tau_L q1_on_dt_R ...
                        tau_R t_P_interp_R t_P_interp_L    
                      
                    
  
elseif (all_IRF == 2)   % to scale accumulation perturbations   
   
  save all_impulse_response_function_values_vol.mat ...
        C_0_all D_0_all index_ldm_L index_ldm_R C_0_ldm_L D_0_ldm_L ...
        C_0_ldm_R D_0_ldm_R test_stop_L test_stop_R time_stop_L time_stop_R ...
        h1_ldm_all_L h1_ldm_L h1_ldm_normalized_L h1_ldm_all_R h1_ldm_R ...
        h1_ldm_normalized_R q1_ldm_L q1_ldm_R q1_ldm_normalized_L ...
        q1_ldm_normalized_R tau_L tau_R t_grid_imp_L t_grid_imp_R ...
        tau_volume_L tau_volume_R alpha1_ldm_L alpha1_ldm_R  
    
    
  % this gets loaded later in convolution.m           
  save q1_on_dt_vol.mat q1_on_dt_L h1_on_dt_L h1_on_dt_R ...
                        tau_L q1_on_dt_R tau_R t_P_interp_R t_P_interp_L
  
end


  
% check to make sure IRF is reasonable:
  if ( (abs(q1_on_dt_L(2)-q1_on_dt_L(1)) < 1e-6) || (abs(q1_on_dt_R(2)-q1_on_dt_R(1)) < 1e-6) )
      disp('Response function incorrect. CHECK response_function_calc.m')
      stop;
  end
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETAILED COMMENTS:

% REFERENCES:
% Nye (1960), Proc Royal Soc, Ser. A, 256 (1287), 559-584
% Nye (1963a), Geophys Jour Royal Astron Soc, 7 (4), 431-456
% Nye (1963b), Proc Royal Soc, Ser. A, 275 (1360), 87-112
% Nye (1965), J Glac, 5(41), 589-606
% Nye (1965), J Glac, 5(41), 567-587 
% Johannesson, et al. (1989), J Glac, 35(121), 355-369


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
