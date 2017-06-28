function [ recalc_impulse_response ] = ...
                         check_recalc_impulse_response( all_IRF, ...
                                                        calculating_jacobians, ...
                                                        force_calc_IRF )


                     
%--------------------------------------------------------------------------
%
% Michelle Koutnik 
% last updated: APRIL 2010, Copenhagen


%--------------------------------------------------------------------------
                     
                     
global solver_method

global itnum_thermomechanical
global pert_distribution
global a_pert




if ((solver_method ~= 4) && (itnum_thermomechanical == 1) && ...
    (calculating_jacobians ~= 1 ) || (force_calc_IRF == 1) )  
% pa engelsk: if *not* using the explicit scheme, if it *is* the first iteration
% of the thermomechanical calculation, and if *not* calculating numerical
% jacobians, then go ahead and...
    


  recalc_impulse_response = 1;    % 1 = calc new impulse response function (IRF)
                                   % 0 = use saved IRF

                                   

  
 % ---------------- 
  if (all_IRF == 1)      % IRF to scale volume change from ice-divide migration
 % ---------------- 
 
 
 % delta-function perturbation:
   pert_distribution       = 'nonuniform';     % DEFAULT
   a_pert                  = 1/10;             %   


   
 
 % --------------------  
  elseif (all_IRF == 2)    % IRF to scale volume change from accumulation changes
 % -------------------- 
  
% (these other options were used for testing)

 
% bdot perturbation over the full domain:
% ---------------------------------------  
   pert_distribution       = 'full_uniform';    % DEFAULT
   
 %  pert_distribution       = 'full_nonuniform';
 %  a_pert                  = 1/2;       % must declare if using nonuniform.
 
 %  pert_distribution       = 'full_linear';
 %  pert_distribution       = 'full';   % interpolate from bdot pattern. 
 
 
% bdot perturbation over only the limited domain:
% -----------------------------------------------
  %  pert_distribution       = 'uniform';  
  
  %  pert_distribution       = 'nonuniform';
  %  a_pert                  = 1/2;        % must declare if using nonuniform.

  %  pert_distribution       = 'linear';

   

   
   
  end   % end if statement on all_IRF
  

  
  
else
    
    
 recalc_impulse_response = 0;      % use saved impulse response function
 
 
 
end   % end if statement on recalc_impulse_response




if (force_calc_IRF == 0)
    
    recalc_impulse_response = 0;

end





% _________________________________________________________________________
% % Keep bits around just in case.

% if ( (i_time == 1) || ((i_time/reset_time) > 1) && (itnum_dT == 1) )
% 
%   if (i_time <= reset_time)
%     i_time_use = i_time;                                % set to time = 1
%   elseif (i_time > reset_time)
%     i_time_use = (reset_time*floor(i_time/reset_time)); % update
%   end 



% If recalculating the impulse response function, pass in a new surface --
% is it ok if it is not in steady state??
