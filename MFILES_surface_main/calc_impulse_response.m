function [ q1_on_dt_L_vol, ...
           q1_on_dt_R_vol, ...
           q1_on_dt_L_div, ...
           q1_on_dt_R_div ] = calc_impulse_response ( h_P, x_P, t_P, ...
                                                      b_dot_P, ...
                                                      x_e, x_w, dx_P, ...
                                                      dx_e, dx_w, ...
                                                      B_P, W_P, slip_P, ...
                                                      A_eff_edges, ...
                                                      Q_out_L_SS, Q_out_R_SS, ...
                                                      Q_external_L, ...
                                                      Q_external_R, ...
                                                      divide_pos_P_SS, ...
                                                      divide_pos_real_SS, ...
                                                      divide_position, ...
                                                      i_time, dt_P, ...
                                                      u_bar_edges, ...
                                                      b_dot_edges_all_t, ...
                                                      b_dot_P_all_t )

                          
% ----------------------------------------------------------------------- %
% Michelle Koutnik
% last updated: March 2010, Copenhagen


% Steps:
% 1) embed the limited-domain in a full model (ldm_extension.m)
%    with a wedge terminus (wedge_terminus.m)
% 2) check full surface holds steady state with wedge terminus
% 3) calculate ice-thickness evolution to impulsive perturbation
% 4) calculate flux response functions

%
% ----------------------------------------------------------------------- %


global dt_impulse
global use_vialov_ext paterson_c_over_a
global ss_return_value
global itnum_max
global pert_distribution
global a_pert
global all_IRF
global perturbation_vol





if (all_IRF == 1)

    
% -------------------------------------------------------------------------
% 1. embed the limited surface in a full domain
% -------------------------------------------------------------------------


 [ x_ldm_P, h_ldm_P, x_div, ...
   x_div_full, index_div_full, ...
   full_x_P, full_x_w, full_x_e, ...
   full_dx_P, full_dx_w, full_dx_e, ...
   full_h_P, full_h_w, full_h_e, ...
   full_S_P, full_S_w, full_S_e, ...
   full_B_P, full_B_w, full_B_e, ...
   full_W_P, full_W_w, full_W_e, ...
   full_slip_P, full_slip_w, full_slip_e, ...
   full_b_dot_P, full_b_dot_w, full_b_dot_e, ...
   full_A_eff_edges, ...
   wedge_length_L, ...
   wedge_length_R, ...
   Q_in_L, Q_in_R, ...
   ablation_rate_L, ...
   ablation_rate_R ] = ldm_extension ( h_P, x_P, ...
                                       x_e, x_w, dx_P, dx_e, dx_w, ...
                                       B_P, W_P, slip_P, ...
                                       A_eff_edges, ...
                                       b_dot_P, divide_position );       
                          
                                   
  save ldm_extension_values.mat  ...
       ...
        x_ldm_P h_ldm_P index_div_full x_div full_x_P full_x_w full_x_e ...
        full_dx_P full_dx_w full_dx_e full_h_P full_h_w full_h_e ...
        full_B_P full_B_w full_B_e ...
        full_W_P full_W_w full_W_e full_S_P full_S_w full_S_e ...
        full_slip_P full_slip_w full_slip_e full_b_dot_P ...
        full_b_dot_w full_b_dot_e full_A_eff_edges ...
        wedge_length_L wedge_length_R ablation_rate_L ablation_rate_R ...
        Q_in_L Q_in_R use_vialov_ext paterson_c_over_a ...
        h_P x_P b_dot_P divide_position
                             
                             
  
    

% Transient calculations with full domain
% ========================================

  full_domain = 1;   % 1 = it is a full domain; 0 = limited domain
  
  
 
  
  
  
% -------------------------------------------------------------------------
% 2. make sure holds steady state with wedge terminus
% -------------------------------------------------------------------------

  disp('Checking that full domain with wedges will hold steady state...')

  
  
  N_time_impulse_ss    = itnum_max;   % maximum number of iterations.
  time_vec_ss          = [ 1: dt_impulse: itnum_max*dt_impulse ];  % dt_impulse is global

  
                                                                
 
% setup accumulation rate, Q_in, and A_eff over time
% ==================================================
  
  b_dot_P_t            = repmat(full_b_dot_P, N_time_impulse_ss, 1); 
  A_eff_edges_xt       = repmat(full_A_eff_edges, N_time_impulse_ss, 1);
  b_dot_edges_t        = repmat([full_b_dot_w full_b_dot_e(end)], ...
                                N_time_impulse_ss, 1);
  b_dot_edges_SS       = b_dot_edges_t(1,:);
   
                             
% initialize values to be filled
% ===============================  
  h_P_xt_ss            = NaN( N_time_impulse_ss, length(full_x_P) );
  h_w_xt_ss            = NaN( N_time_impulse_ss, length(full_x_w) );
  h_e_xt_ss            = NaN( N_time_impulse_ss, length(full_x_e) );
  S_P_xt_ss            = NaN( N_time_impulse_ss, length(full_x_P) );
  S_w_xt_ss            = NaN( N_time_impulse_ss, length(full_x_w) );
  S_e_xt_ss            = NaN( N_time_impulse_ss, length(full_x_e) );
  wedge_length_L_ss    = NaN( N_time_impulse_ss, 1 );
  wedge_length_R_ss    = NaN( N_time_impulse_ss, 1 );
  %flux_edges_kin_xt_ss = NaN( N_time_impulse_ss, length(full_x_P)+1 );
  flux_edges_dyn_xt_ss = NaN( N_time_impulse_ss, length(full_x_P)+1 );
  divide_pos_P_ss      = NaN( N_time_impulse_ss, 1 );
  divide_pos_real_ss   = NaN( N_time_impulse_ss, 1 );
  
  h_P_xt_ss(1,:)       = full_h_P;    % fill first timestep with known value
  h_w_xt_ss(1,:)       = full_h_w;
  h_e_xt_ss(1,:)       = full_h_e;
  S_P_xt_ss(1,:)       = full_S_P;    % fill first timestep with known value
  S_w_xt_ss(1,:)       = full_S_w;
  S_e_xt_ss(1,:)       = full_S_e;
  wedge_length_L_ss(1) = wedge_length_L;
  wedge_length_R_ss(1) = wedge_length_R;
  Q_out_L              = 0;            % by definition for full model
  Q_out_R              = 0;
  divide_pos_P_ss(1)   = x_div_full;   % the new divide is at x=0, but don't need divide position here.
  divide_pos_real_ss(1) = x_div_full;
  
  
  
% TIME loop
% =========
       ii = 1;   % start at time=2
       calc_ss_return_value = 9999;
       disp(['Steady-state return value is ' num2str(ss_return_value)])
  
       
   while ((calc_ss_return_value > ss_return_value) && (ii < itnum_max))
       
       ii = ii+1;   % start at time = 2
          
       
 [ h_P_t, h_w_t, h_e_t, ...
   S_P_t, S_w_t, S_e_t, ...
   wedge_length_L, wedge_length_R, ...
   flux_edges_dyn_t, ...
   vol_change_L_t, vol_change_R_t, ...
   div_change_t, ...
   Q_out_L_t, ...
   Q_out_R_t, ...
   divide_pos_t, ...
   divide_real_t, ...
   Q_external_check_L, ...
   Q_external_check_R] =  implicit_solver2( h_P_xt_ss(ii-1,:), S_P_xt_ss(ii-1,:), ...
                                        divide_pos_P_ss(ii-1), ...
                                        divide_pos_real_ss(ii-1), ...
                                        full_x_P, full_x_w, full_x_e, ...
                                        full_dx_P, full_dx_w, full_dx_e, ...
                                        full_B_P, full_B_w, full_B_e, ...
                                        full_W_P, full_W_w, full_W_e, ...
                                        full_slip_P, full_slip_w, full_slip_e, ...
                                        A_eff_edges_xt(ii-1:ii,:), ...
                                        b_dot_P_t(ii-1:ii,:), ...
                                        b_dot_edges_SS, ...
                                        b_dot_edges_t(ii-1:ii,:), ...
                                        divide_pos_P_SS, ...
                                        divide_pos_real_SS, ...
                                        Q_out_L_SS, Q_out_R_SS, ...
                                        Q_external_L, Q_external_R, ...
                                        Q_out_L, Q_out_R, full_domain, ...
                                        t_P, dt_impulse, i_time, 0,0,0,0 );                             
                                                                 % send zeros as placeholders.
  
  h_P_xt_ss(ii,:)            = h_P_t;
  h_w_xt_ss(ii,:)            = h_w_t;
  h_e_xt_ss(ii,:)            = h_e_t;
  S_P_xt_ss(ii,:)            = S_P_t;
  S_w_xt_ss(ii,:)            = S_w_t;
  S_e_xt_ss(ii,:)            = S_e_t;
  wedge_length_L_ss(ii)      = wedge_length_L;
  wedge_length_R_ss(ii)      = wedge_length_R;
%  flux_edges_kin_xt_ss(ii,:) = flux_edges_kin_t;
  flux_edges_dyn_xt_ss(ii,:) = flux_edges_dyn_t;
  divide_pos_P_ss(ii)        = divide_pos_t;
  divide_pos_real_ss(ii)     = divide_real_t;
  
  % vol_change_L_t, vol_change_R_t, div_change_t, Q_out_L_t, Q_out_R_t
  % are not generated in calculations for full-domain model.
  
  
  
 % keep running until difference in surface is small
 % (until it holds steady state)
 % =================================================
  calc_ss_return_value = max(abs(h_P_xt_ss(ii,:) - h_P_xt_ss(ii-1,:)))

  
  
  end   % while loop                

  
  
  if (ii >= itnum_max)
      disp('BOMBED OUT! In calc_impulse_response.m -- Full surface with wedge did not hold steady state')
      save ss_values_attempt.mat
      stop;
  else
      disp(['Surface with wedge holds steady state. Required ' num2str(ii) ' iterations.'])
      index_end = ii;
  end
      
      
  
   save ss_values.mat ...
        ...
          h_P_xt_ss S_P_xt_ss time_vec_ss dt_impulse wedge_length_L_ss ...
          h_w_xt_ss h_e_xt_ss S_w_xt_ss S_e_xt_ss ...
          full_x_P full_x_w full_x_e ...
          flux_edges_dyn_xt_ss ...
          b_dot_P_t N_time_impulse_ss b_dot_edges_SS ...
          index_end wedge_length_R_ss divide_pos_P_ss ...
          divide_pos_real_ss full_domain
                
 



                  
elseif (all_IRF > 1)  % the saved files can be used.               
                  
                  
 load ldm_extension_values.mat
 load ss_values.mat
 full_domain = 1;


end    % if statement on all_IRF





% -------------------------------------------------------------------------
% 3. calculate surface evolution to impulsive perturbation
% -------------------------------------------------------------------------


disp( ' ')
disp('Calculating evolution of extended domain with wedges to impulsive perturbation...')
disp(['The distribution of the accumulation perturbation is ' pert_distribution '.'])

  
  N_time_impulse = itnum_max;   % maximum number of iterations.
  time_vec       = [ 1: dt_impulse: itnum_max*dt_impulse ];
  % uses global value 'dt_impulse' here!
  
  
  time_perturb   = 2;           % add perturbation at second timestep
 
 
  full_x_edges      = [full_x_w(1) full_x_e];
  full_x_edges_norm = [full_x_w(1) full_x_e] / full_x_w(1);
  full_x_P_norm     = full_x_P / full_x_P(1);
  
  
% pull out the index of the limited-domain boundary
% on the full-domain grid
% =======================
  index_ldm_L = find((full_x_P <= (x_P(1) - x_div)), 1, 'last');
  index_ldm_R = find((full_x_P <= (x_P(end) - x_div)), 1, 'last'); 
  
  
% initialize perturbation vectors
% ================================
  perturbation_P     = zeros(size(full_x_P));
  perturbation_edges = zeros(size(full_x_edges));
  
  

 % extra points (ep) beyond LD
 % ============================
   ep = 2;
   
 
   
  
% lots of options here, but in standard runs the only ones used are
% "full_uniform" (to scale volume changes due to accumulation changes)
% and "nonuniform" (to scale volume changes due to divide migration)
   
   
% uniform perturbation:
% ---------------------  
  if (strcmp(pert_distribution, 'uniform') == 1)
      
     perturbation_edges( index_ldm_L-ep:index_ldm_R+ep ) = (perturbation_vol / ...
                        length(full_x_edges(index_ldm_L-ep:index_ldm_R+ep)));
                         
                         
     perturbation_P( index_ldm_L-ep:index_ldm_R+ep ) = (perturbation_vol / ...
                              length(full_x_P(index_ldm_L-ep:index_ldm_R+ep)));

                          
 
% full uniform perturbation:
% --------------------------
  elseif (strcmp(pert_distribution, 'full_uniform') == 1)
      
     perturbation_edges = ones(size(full_x_edges)) * (perturbation_vol / length(full_x_edges));
                         
     perturbation_P = ones(size(full_x_P)) * (perturbation_vol / length(full_x_P));                          
     
     

%  nonuniform perturbation:
%  ------------------------     
  elseif (strcmp(pert_distribution, 'nonuniform') == 1)                         

   pert_edges_temp    = (1/(a_pert*sqrt(pi))) * ...
                        exp(-(full_x_edges_norm(index_ldm_L-ep:index_ldm_R+ep)).^2 / (a_pert^2));
   pert_edges_norm    = pert_edges_temp / sum(pert_edges_temp);
   
   perturbation_edges(index_ldm_L-ep:index_ldm_R+ep) = ...
                        pert_edges_norm * perturbation_vol;
  
   pert_P_temp        = (1/(a_pert*sqrt(pi))) * ...
                        exp(-(full_x_P_norm(index_ldm_L-ep:index_ldm_R+ep)).^2 / (a_pert^2));
   pert_P_norm        = pert_P_temp / sum(pert_P_temp);
   perturbation_P(index_ldm_L-ep:index_ldm_R+ep) = ...
                        pert_P_norm * perturbation_vol;     

                    
                    
%  full nonuniform perturbation:
%  -----------------------------     
  elseif (strcmp(pert_distribution, 'full_nonuniform') == 1)                         

   pert_edges_temp    = (1/(a_pert*sqrt(pi))) * exp(-(full_x_edges_norm).^2 / (a_pert^2));
   pert_edges_norm    = pert_edges_temp / sum(pert_edges_temp);
   perturbation_edges = pert_edges_norm * perturbation_vol;
  
   pert_P_temp        = (1/(a_pert*sqrt(pi))) * exp(-(full_x_P_norm).^2 / (a_pert^2));
   pert_P_norm        = pert_P_temp / sum(pert_P_temp);
   perturbation_P     = pert_P_norm * perturbation_vol;       
                    
                    
   
%  linearly varing:
%  ----------------   
  elseif (strcmp(pert_distribution, 'linear') == 1)                   
    
    edge_diff = (perturbation_vol / ...
                  length(full_x_edges(index_ldm_L-ep:index_ldm_R+ep)));
    P_diff    = (perturbation_vol / ...
                  length(full_x_P(index_ldm_L-ep:index_ldm_R+ep)));
    
    % want the perturbation to be positive
    perturbation_edges(index_ldm_L-ep:index_ldm_R+ep)  = ...
                          ((edge_diff) + ...
                          [-1:(2/(length(full_x_edges(index_ldm_L-ep:index_ldm_R+ep))-1)):1] * edge_diff);
                      
    perturbation_P(index_ldm_L-ep:index_ldm_R+ep) = ...
                          (P_diff) + ...
                          ([-1: (2/(length(full_x_P(index_ldm_L-ep:index_ldm_R+ep))-1)): 1] * P_diff);
  

%  full linearly varing:
%  ---------------------   
  elseif (strcmp(pert_distribution, 'full_linear') == 1)                   
    
    edge_diff = (perturbation_vol / length(full_x_edges));
    P_diff    = (perturbation_vol / length(full_x_P));
    
    % want the perturbation to be positive
    perturbation_edges  = ((edge_diff) + ...
                          [-1:(2/(length(full_x_edges)-1)):1] * edge_diff);
                      
    perturbation_P = (P_diff) + ...
                     ([-1: (2/(length(full_x_P)-1)): 1] * P_diff);
  
                      
                      
% over full domain:
% -----------------
  elseif (strcmp(pert_distribution, 'full') == 1)                   
                   
% for the given mass-balance pattern, interpolate as done in test with full model
% get the pattern from bdot across limited domain

   t_max_bdot = find(t_P(max(b_dot_edges_all_t(:,1))));

   perturbation_edges = interp1( [x_w(1) x_e]-x_div, ...
                             (b_dot_edges_all_t(t_max_bdot,:)-b_dot_edges_all_t(1,:)), ...
                              [full_x_w(1) full_x_e], 'linear', 'extrap');   
                     
% make sure it sums to perturbation volume:
   addfactor = (perturbation_vol - ...
                (sum(perturbation_edges)))/length(perturbation_edges);
   perturbation_edges = perturbation_edges + addfactor;
             
   
   perturbation_P  = interp1( x_P-x_div, (b_dot_P_all_t(t_max_bdot,:)-b_dot_P_all_t(1,:)), ...
                              full_x_P, 'linear', 'extrap' );
 
   addfactor = (perturbation_vol - ...
                (sum(perturbation_P)))/length(perturbation_P);
   perturbation_P = perturbation_P + addfactor;
      
   
  end  % if on pert_distribution  
  




                    
  

 % setup accumulation rate, Q_in, and A_eff over time
 % ==================================================
 
  b_dot_P_t      = repmat(full_b_dot_P, N_time_impulse, 1); 
  b_dot_edges_t  = repmat([full_b_dot_w full_b_dot_e(end)], ...
                           N_time_impulse, 1);
  A_eff_edges_xt = repmat(full_A_eff_edges, N_time_impulse, 1);
  
  
  
  % now, hit it with accumulation perturbation! bam.
  b_dot_P_t(time_perturb,:)     = b_dot_P_t(time_perturb,:) + perturbation_P;   
  b_dot_edges_t(time_perturb,:) = b_dot_edges_t(time_perturb,:) + perturbation_edges;
 
  
                             
                             
  % initialize values to be filled
  % ===============================
  h_P_xt                 = NaN( N_time_impulse, length(full_x_P) );
  h_w_xt                 = NaN( N_time_impulse, length(full_x_w) );
  h_e_xt                 = NaN( N_time_impulse, length(full_x_e) );
  S_P_xt                 = NaN( N_time_impulse, length(full_x_P) );
  S_w_xt                 = NaN( N_time_impulse, length(full_x_w) );
  S_e_xt                 = NaN( N_time_impulse, length(full_x_e) );
  wedge_length_L_t       = NaN( N_time_impulse, 1 );
  wedge_length_R_t       = NaN( N_time_impulse, 1 );
  flux_edges_kin_xt      = NaN( N_time_impulse, length(full_x_P)+1 );
  flux_edges_dyn_xt      = NaN( N_time_impulse, length(full_x_P)+1 );
  divide_pos_P_t         = NaN( N_time_impulse, 1 );
  divide_pos_real_t      = NaN( N_time_impulse, 1 );
  
  h_P_xt(1,:)            = h_P_xt_ss(index_end,:);   % fill first timestep with SS value
  h_w_xt(1,:)            = h_w_xt_ss(index_end,:);
  h_e_xt(1,:)            = h_e_xt_ss(index_end,:);
  S_P_xt(1,:)            = S_P_xt_ss(index_end,:);   % fill first timestep with SS value
  S_w_xt(1,:)            = S_w_xt_ss(index_end,:);
  S_e_xt(1,:)            = S_e_xt_ss(index_end,:);
  wedge_length_L_t(1,:)  = wedge_length_L_ss(index_end);
  wedge_length_R_t(1,:)  = wedge_length_R_ss(index_end);
%  flux_edges_kin_xt(1,:) = flux_edges_kin_xt_ss(index_end,:);
  flux_edges_dyn_xt(1,:) = flux_edges_dyn_xt_ss(index_end,:);
  divide_pos_P_t(1,:)    = divide_pos_P_ss(index_end);
  divide_pos_real_t(1,:) = divide_pos_real_ss(index_end);
  Q_out_L                = 0;          % by definition for full model
  Q_out_R                = 0;
  
  
% TIME loop
% =========
 
       ii = 1;   % start at time=2
       calc_ss_return_value = 9999;

     
       disp(['The return value is ' num2str(ss_return_value)])
  
  
   while ((calc_ss_return_value > ss_return_value) && (ii < itnum_max))
       
       ii = ii+1;
       
     
 [ h_P_t, h_w_t, h_e_t, ...
   S_P_t, S_w_t, S_e_t, ...
   wedge_length_L, wedge_length_R, ...
   flux_edges_dyn_t, ...
   vol_change_L_t, vol_change_R_t, ...
   div_change_t, ...
   Q_out_L_t, ...
   Q_out_R_t, ...
   divide_pos_t, ...
   divide_real_t ] =  implicit_solver2( h_P_xt(ii-1,:), ...
                                        S_P_xt(ii-1,:), divide_pos_P_t(ii-1), ...
                                        divide_pos_real_t(ii-1), ...
                                        full_x_P, full_x_w, full_x_e, ...
                                        full_dx_P, full_dx_w, full_dx_e, ...
                                        full_B_P, full_B_w, full_B_e, ...
                                        full_W_P, full_W_w, full_W_e, ...
                                        full_slip_P, full_slip_w, full_slip_e, ...
                                        A_eff_edges_xt(ii-1:ii,:), ...
                                        b_dot_P_t(ii-1:ii,:), ...
                                        b_dot_edges_SS, b_dot_edges_t(ii-1:ii,:), ...
                                        divide_pos_P_SS, ...
                                        divide_pos_real_SS, ...
                                        Q_out_L_SS, Q_out_R_SS, ...
                                        Q_external_L, Q_external_R, ...
                                        Q_out_L, Q_out_R, ...
                                        full_domain, t_P, ...
                                        dt_impulse, i_time, 0,0,0,0 );
                                % send zeros for the response function
                                

  h_P_xt(ii,:)            = h_P_t;
  h_w_xt(ii,:)            = h_w_t;
  h_e_xt(ii,:)            = h_e_t;
  S_P_xt(ii,:)            = S_P_t;
  S_w_xt(ii,:)            = S_w_t;
  S_e_xt(ii,:)            = S_e_t;
  wedge_length_L_t(ii)    = wedge_length_L;
  wedge_length_R_t(ii)    = wedge_length_R;
 % flux_edges_kin_xt(ii,:) = flux_edges_kin_t;
  flux_edges_dyn_xt(ii,:) = flux_edges_dyn_t;
  divide_pos_P_t(ii)      = divide_pos_t;
  divide_pos_real_t(ii)   = divide_real_t;
  
  
  
 % keep running until difference in surface is small
 % (until returns to steady state)
 % =================================================
   calc_ss_return_value = max(abs(h_P_xt(ii,:) - h_P_xt(ii-1,:)));  % accounts for both sides
  
  end   % while loop                
  
  
  if (ii >= itnum_max)
      disp('BOMBED OUT! In calc_impulse_response.m -- Full surface with wedge did not find IRF')
      save impulse_values_attempt.mat
      stop;
  else
      disp(['Surface returned to steady state after perturbation. Required ' num2str(ii) ' iterations.'])
  end
      
                  
  end_value         = ii;    % keep only values that are not NaN
  
  time_vec          = time_vec(1:end_value);
  h_P_xt            = h_P_xt(1:end_value,:);
  h_e_xt            = h_e_xt(1:end_value,:);
  h_w_xt            = h_w_xt(1:end_value,:);
 % flux_edges_kin_xt = flux_edges_kin_xt(1:end_value,:);
  flux_edges_dyn_xt = flux_edges_dyn_xt(1:end_value,:);
  b_dot_P_t         = b_dot_P_t(1:end_value,:);
  
  
  save impulse_values.mat ...
       ...
          full_x_P full_x_w h_P_xt full_x_e h_e_xt h_w_xt...
          S_P_xt S_w_xt S_e_xt time_vec dt_impulse wedge_length_L_t ...
          time_perturb flux_edges_dyn_xt b_dot_P_t ...
          N_time_impulse perturbation_vol perturbation_P ...
          perturbation_edges end_value wedge_length_R_t ...
          divide_pos_P_t divide_pos_real_t index_ldm_L index_ldm_R

                      






% -------------------------------------------------------------------------
% 4. calculate flux response function
% -------------------------------------------------------------------------


%  response_function_calc_dynflux;   % calculates q1 from dynamic flux using
                                    % numerical values from implicit solver

 response_function_calc_q1h1;      % calculates q1 from h1 using kinematic
                                    % wave theory

  
  
  
if (all_IRF == 2)    % both response functions have been calculated.

   load q1_on_dt_vol.mat

   q1_on_dt_L_vol = q1_on_dt_L;
   q1_on_dt_R_vol = q1_on_dt_R;

   load q1_on_dt_div.mat

   q1_on_dt_L_div = q1_on_dt_L;
   q1_on_dt_R_div = q1_on_dt_R;


else   % must pass out some values. 
    
    q1_on_dt_L_vol = NaN;
    q1_on_dt_R_vol = NaN;
    q1_on_dt_L_div = NaN;
    q1_on_dt_R_div = NaN;
    
end






