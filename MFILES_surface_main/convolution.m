function [ Q_out_L, ...
           Q_out_R, ...
           Q_external_check_L, ...
           Q_external_check_R ] = convolution ( vol_change_L, vol_change_R, ...
                                     div_change, divide_pos_real_SS, ...
                                     x_divide_real, Q_out_L_SS, Q_out_R_SS, ...
                                     Q_external_L, Q_external_R, ...
                                     i_time, t_P, dt_P, ...
                                     q1_on_dt_L_vol, q1_on_dt_R_vol, ...
                                     q1_on_dt_L_div, q1_on_dt_R_div )


% -------------------------------------------------------------------------
% Michelle Koutnik
% last updated: April 2008


% want the integral from 0 to tau of q1 = 1
% this multiplies the volume perturbation

% get the time series of spatially averaged volume input perturbations
%  Q(x_flux(j), t_flux ) - Q(x_flux(j), t_flux-1) = 
%       Q_in( t_flux ) - Q_in( t_flux-1 )
%       + integral_{x_div(t)}^{x_flux(j)}[ b_dot(x,t) - b_dot(x,t-1)] W(x)dx
 
% -------------------------------------------------------------------------


global N_t_mesh
global conv_prod_matrix_vol_L    % try this as a global variable.
global conv_prod_matrix_div_L
global conv_prod_matrix_vol_R
global conv_prod_matrix_div_R



if (i_time == 1)
  % I don't think i_time = 1 is ever calculated, but just in case...
  
  Q_out_L = Q_out_L_SS;
  Q_out_R = Q_out_R_SS;
  
  
  
elseif (i_time >= 2)    
    
   
 % load q1_on_dt_vol.mat
  
  
  q1_vol_use_L = q1_on_dt_L_vol * vol_change_L;   % this is for one timestep
  q1_vol_use_R = q1_on_dt_R_vol * vol_change_R;   % this is for one timestep

  

 % load q1_on_dt_div.mat 
   
 
 
  if (divide_pos_real_SS < x_divide_real)      % divide moves to the right.

      q1_div_use_L = q1_on_dt_L_div * (-div_change);  % increases
      q1_div_use_R = q1_on_dt_R_div * (-div_change);  % decreases
      
  elseif (divide_pos_real_SS > x_divide_real)  % divide moves to the left

      q1_div_use_L = q1_on_dt_L_div * (div_change);    % decreases
      q1_div_use_R = q1_on_dt_R_div * (div_change);    % increases
      
  elseif (divide_pos_real_SS == x_divide_real)
      q1_div_use_L = 0;
      q1_div_use_R = 0;
  
  end

  
  
% ------------
% convolution
% ------------

% tau is the time when the impulse response function has returned to
% steady-state value.  only add contributions over tau time.

if (i_time == 2)  
    
% fill matrix of values for every timestep
% ========================================
   conv_prod_matrix_vol_L = NaN( N_t_mesh, N_t_mesh );
   conv_prod_matrix_div_L = NaN( N_t_mesh, N_t_mesh );
   conv_prod_matrix_vol_R = NaN( N_t_mesh, N_t_mesh );
   conv_prod_matrix_div_R = NaN( N_t_mesh, N_t_mesh );
     
   
 % there is no volume change for the first timestep, 
 % so fill the first row with zeros
 % ================================
   conv_prod_matrix_vol_L(1, 1:length(q1_on_dt_L_vol)) = zeros(1,length(q1_on_dt_L_vol));
   conv_prod_matrix_div_L(1, 1:length(q1_on_dt_L_div)) = zeros(1,length(q1_on_dt_L_div));
   conv_prod_matrix_vol_R(1, 1:length(q1_on_dt_R_vol)) = zeros(1,length(q1_on_dt_R_vol));
   conv_prod_matrix_div_R(1, 1:length(q1_on_dt_R_div)) = zeros(1,length(q1_on_dt_R_div));

end  % on if i_time == 2
      
      
   


% fill for the current timestep.
% ==============================
   conv_prod_matrix_vol_L(i_time,1:length(q1_on_dt_L_vol)) = q1_vol_use_L;
   conv_prod_matrix_div_L(i_time,1:length(q1_on_dt_L_div)) = q1_div_use_L;
   conv_prod_matrix_vol_R(i_time,1:length(q1_on_dt_R_vol)) = q1_vol_use_R;
   conv_prod_matrix_div_R(i_time,1:length(q1_on_dt_R_div)) = q1_div_use_R;

  


% want all the values in each column over tau time
% e.g. time=1, (1,1); time=2, (1,2)+(2,1); time=3, (1,3)+(2,2)+(3,1)
% ==================================================================


% for the left:
% =============
% if ( (i_time*dt_P) <= tau_L )
%   start_time = 1;
% elseif ( (i_time*dt_P) > tau_L )
%   start_time = 1+(i_time+1-round(tau_L/dt_P));
% end

start_time = 1;   % filled the rest with zeros.


% initialize values and then sum in loop:
    Q_add_vol_L = 0;
    Q_add_div_L = 0;
     
   for i = i_time:-1:start_time  
        Q_add_vol_L = Q_add_vol_L + conv_prod_matrix_vol_L(i, i_time+1-i);
        Q_add_div_L = Q_add_div_L + conv_prod_matrix_div_L(i, i_time+1-i);
   end
   
   
   
% for the right:
% ==============
% if ( (i_time*dt_P) <= tau_R )
%   start_time = 1;
% elseif ( (i_time*dt_P) > tau_R )
%   start_time = 1+(i_time+1-round(tau_R/dt_P));
% end


% initialize values and then sum in loop:
    Q_add_vol_R = 0;
    Q_add_div_R = 0;
    
   for i = i_time:-1:start_time  % NEED TO CHANGE IN CASE size(q1_on_dt) < N_t_mesh
        Q_add_vol_R = Q_add_vol_R + conv_prod_matrix_vol_R(i, i_time+1-i);
        Q_add_div_R = Q_add_div_R + conv_prod_matrix_div_R(i, i_time+1-i);
   end   
   
   

% add convolution product to steady-state fluxes
% ==============================================

% % Attempt to balance changes in ice thickness with changes in Qext, 8 August 2014
% % disp('Forcing value of Qext in convolution.m!!')
% Q_external_check_L = -Q_add_vol_L - Q_add_div_L;
% Q_external_check_R = -Q_add_vol_R - Q_add_div_R;
% Q_external_L = Q_external_check_L;
% Q_external_R = Q_external_check_R;


% Don't change the value:
 Q_external_check_L = Q_external_L;
 Q_external_check_R = Q_external_R;


  Q_out_L = Q_out_L_SS + Q_external_L + Q_add_vol_L + Q_add_div_L;
  
  Q_out_R = Q_out_R_SS + Q_external_R + Q_add_vol_R + Q_add_div_R;


   
  
  
% Make sure NaN didn't stick around...  
if( (isnan(Q_out_L) == 1) || (isnan(Q_out_R) == 1) )
    disp('NaN in convolution.m; Check q1_on_dt.mat from response_function_calc.m')
    stop;
end



end   % on i_time >= 2



