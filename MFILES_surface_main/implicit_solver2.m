function [ h_P_t, h_w_t, h_e_t, ...
           S_P_t, S_w_t, S_e_t, ...
           flux_edges_dyn_t     ] = implicit_solver2( h_P, S_P, ...
                                                 x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e, ...
                                                 B_P, B_w, B_e, ...
                                                 W_P, W_w, W_e,  ...
                                                 E_P, E_w, E_e, ...
                                                 fs_P, fs_w, fs_e, ...
                                                 A_eff_edges, b_dot_P, ...
                                                 b_dot_edges_SS, b_dot_edges, ...
                                                 Q_out_L_SS, Q_out_R_SS, ...
                                                 Q_external_L, Q_external_R, ...
                                                 S_at_GL, t_P, ...
                                                 dt_P, i_time, ...
                                                 deformation_only, deformation_plus_sliding, sliding_only)

                                                          
% -------------------------------------------------------------------------
% Michelle Koutnik
% last updated: March 2010, Copenhagen

% Find surface evolution using implicit method (#2)
% This solver accounts for "diffusion" only.

% can work for full model, or calculates boundary-flux values using
% impulse-response functions for the limited-domain model.



% NOTES:
% - underrelaxation must be used (beta2=0.1, set in global_variables.m)

% -------------------------------------------------------------------------


global rho_ice  g  n
global itnum_max 
global res_stop2
global calculate_residual

% global deformation_only deformation_plus_sliding sliding_only ...
%        deformation_sliding_lateraldrag deformation_sliding_longstress ...
%        deformation_sliding_lateraldrag_longstress
   

% Initialize values
% =================
  flow_constant   = ((2*A_eff_edges(1,:))/(n+2)) * (rho_ice*g)^n;  
  flow_constant_t = ((2* A_eff_edges(2,:))/(n+2)) * (rho_ice*g)^n;
%  from the shallow ice approximation  


  N_x = length(x_P);   % number of points


  
% known values
% ============
  [ S_w, S_e ] = get_edge_values_quadratic( S_P, x_P, x_w, x_e, ...
                                            dx_P, dx_w, dx_e );
 % S_w(1) = S_at_GL(1);
                                        
  [ h_w, h_e ] = get_edge_values_quadratic( h_P, x_P, x_w, x_e, ...
                                            dx_P, dx_w, dx_e );                                       
    
  [ dS_dx_w, dS_dx_e ] = get_gradient_values( S_P, x_P, dx_P );
  

   
 % make a guess of the unknown value
 % =================================
   S_P_t     = S_P;   
   h_P_t     = h_P;


% -------------------------------------------------
% Iteratively find h_P_t that satisfies continuity
% -------------------------------------------------
    
   max_res = 999999;
   itnum = 0;
   
% iterate until changes are smaller than a threshold value (res_stop2)   
  while ((max_res > res_stop2) && (itnum < itnum_max))
     itnum = itnum+1;
                    
     
  [ S_w_t, S_e_t ] = get_edge_values_quadratic( S_P_t, x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );
  %  S_w_t(1) = S_at_GL(2);                   
                                             
  [ h_w_t, h_e_t ] = get_edge_values_quadratic( h_P_t, x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );                                                                 
  
  [ dS_dx_w_t, dS_dx_e_t ] = get_gradient_values( S_P_t, x_P, dx_P );     
      

  
% Call solver
% =========== 

if (deformation_only == 1)
    
   [ S_P_t_guess, ...
     S_P_t, ...
     flux_edges_dyn_t ] = solver( flow_constant, flow_constant_t, ...
                                  h_P, h_w, h_e, h_P_t, h_w_t, h_e_t, ...
                                  S_P, S_w, S_e, S_P_t, S_w_t, S_e_t, ...
                                  dS_dx_w, dS_dx_e, dS_dx_w_t, dS_dx_e_t, ...
                                  B_P, B_w, B_e, W_P, W_w, W_e, ...
                                  x_P, x_w, x_e, ...
                                  dx_P, dx_w, dx_e, dt_P, ...
                                  b_dot_P, b_dot_edges, ...
                                  S_at_GL, N_x, Q_out_R_SS, ...
                                  Q_external_L, Q_external_R, ...
                                  E_w, E_e, fs_w, fs_e, ...
                                  deformation_only, deformation_plus_sliding, sliding_only);
                              
elseif (deformation_plus_sliding == 1)
              
   [ S_P_t_guess, ...
     S_P_t, ...
     flux_edges_dyn_t ] = solver_deformation_sliding( ...
                                  flow_constant, flow_constant_t, ...
                                  h_P, h_w, h_e, h_P_t, h_w_t, h_e_t, ...
                                  S_P, S_w, S_e, S_P_t, S_w_t, S_e_t, ...
                                  dS_dx_w, dS_dx_e, dS_dx_w_t, dS_dx_e_t, ...
                                  B_P, B_w, B_e, W_P, W_w, W_e, ...
                                  x_P, x_w, x_e, ...
                                  dx_P, dx_w, dx_e, dt_P, ...
                                  b_dot_P, b_dot_edges, ...
                                  S_at_GL, N_x, Q_out_R_SS, ...
                                  Q_external_L, Q_external_R, ...
                                  E_w, E_e, fs_w, fs_e, ...
                                  deformation_only, deformation_plus_sliding, sliding_only);
   % Should dS/dx be negative here??
                              
elseif (sliding_only == 1)                              
        
     [ S_P_t_guess, ...
     S_P_t, ...
     flux_edges_dyn_t ] = solver_sliding_only( ...
                                  flow_constant, flow_constant_t, ...
                                  h_P, h_w, h_e, h_P_t, h_w_t, h_e_t, ...
                                  S_P, S_w, S_e, S_P_t, S_w_t, S_e_t, ...
                                  dS_dx_w, dS_dx_e, dS_dx_w_t, dS_dx_e_t, ...
                                  B_P, B_w, B_e, W_P, W_w, W_e, ...
                                  x_P, x_w, x_e, ...
                                  dx_P, dx_w, dx_e, dt_P, ...
                                  b_dot_P, b_dot_edges, ...
                                  S_at_GL, N_x, Q_out_R_SS, ...
                                  Q_external_L, Q_external_R, ...
                                  E_w, E_e, fs_w, fs_e, ...
                                  deformation_only, deformation_plus_sliding, sliding_only);
    
    
% elseif (deformation_sliding_lateraldrag == 1)
%     
%     
% elseif (deformation_sliding_long_stress == 1)
    

end

% %   Q_out_L_t = flux_edges_dyn_t(1);
% %   Q_out_R_t = flux_edges_dyn_t(end);
  
 
                              
% value of the residual
% =====================       
  max_res = max(abs(S_P_t - S_P_t_guess));

  % use the value of h_P_t returned from solver.
  % stop when difference between iterations becomes small.

  

 end   % while loop on max_res


  
 

% send alert if iterations were cutoff
% ====================================
  if (itnum == itnum_max)
      disp([' ']);
      disp(['MAXIMUM ITERATIONS REACHED!']);
      disp(['maxed out iterating on residual, implicit_solver2']);
      disp(['residual = ' num2str(max_res)]);
  end


  
  
% recalculate one last time with S_P_t
% ====================================
 [ S_w_t, S_e_t ] = get_edge_values_quadratic( S_P_t, x_P, x_w, x_e, ...
                                                dx_P, dx_w, dx_e );
 %  S_w_t(1) = S_at_GL(2);
                                            

  h_P_t = S_P_t - B_P;
  
  
  if (isempty(find(sign(h_P_t) == -1)) ~= 1)
  disp(' thickness less than zero!')
  stop
  end
  
      
  [ h_w_t, h_e_t ] = get_edge_values_quadratic( h_P_t, x_P, x_w, x_e, ...
                                                dx_P, dx_w, dx_e );
  
                                            
 
% update dynamic flux:
% ====================
if (deformation_only == 1)
    
    K_w_t = E_w .* flow_constant_t(1:end-1) .* W_w .* (h_w_t.^(n+2)) ...
            .* ( (dS_dx_w_t.^2).^( (n-1)/2 ) );
    K_e_t = E_e .* flow_constant_t(2:end) .* W_e .* (h_e_t.^(n+2)) ...
            .* ( (dS_dx_e_t.^2).^( (n-1)/2 ) );
        
elseif (sliding_only == 1)
    
    K_w_t = ( ( (fs_w .* ((rho_ice*g)^n .* h_w_t.^(n) .* W_w)) ) ...
            .* ( (dS_dx_w_t.^2).^( (n-1)/2 ) ) );    

    K_e_t = ( ( (fs_e .* ((rho_ice*g)^n .* h_e_t.^(n) .* W_e)) ) ...
            .* ( (dS_dx_e_t.^2).^( (n-1)/2 ) ) );           
    
elseif (deformation_plus_sliding == 1)
    
    K_w_t = ( (( E_w .* flow_constant(1:end-1) .* W_w .* (h_w_t.^(n+2))) ...
                + (fs_w .* ((rho_ice*g)^n .* h_w_t.^(n) .* W_w)) ) ...
            .* ( (dS_dx_w_t.^2).^( (n-1)/2 ) ) );    

    K_e_t = ( (( E_e .* flow_constant(1:end-1) .* W_e .* (h_e_t.^(n+2))) ...
                + (fs_e .* ((rho_ice*g)^n .* h_e_t.^(n) .* W_e)) ) ...
            .* ( (dS_dx_e_t.^2).^( (n-1)/2 ) ) );        
        
end
        
  Q_w_t    = - sign(dS_dx_w_t) .* K_w_t .* abs(dS_dx_w_t);   % flux at western edges at time=time-1
  Q_e_t    = - sign(dS_dx_e_t) .* K_e_t .* abs(dS_dx_e_t);   % flux at eastern edges at time=time-1
 
  
  flux_edges_dyn_t = [Q_w_t(1) Q_e_t];
   

 
 
 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETAILED COMMENTS:

%      1.5-D ice sheet evolution, solved using the 
%      Finite Volume Method with "implicit" time stepping,
%      except the "coefficients" u(x,y) and D(x,y) are time dependent
%      and are evaluated using h and dh/dx from old time step. 
%      Based on Patankar,S. 1980  Numerical Heat Transfer and Fluid Flow
%      Needs m-files:
%            fvm_soln.m               fvm implicit equation solver
%            velo.m                   ice flow velocity
%            b_dot.m                  accumulation rate pattern
%            A.m                      upwinding function
%            F_upwind.m               upwinding function

%  Ice Cap equation
%            dh/dt + div( U_bar h ) = b_dot                          (1)
%
%       h(x,y)     = ice cap thickness
%       b_dot(x,y) = net mass balance
%       U_bar(x,y) = (u_bar,v_bar)  depth-averaged velocity components
%
% Assumed Form of Velocity Field
%     This uses parallel-sided slab depth-dependent shape function phi(z),
%     which is then vertically integrated to give horizontal ice transport.
%     
%        u(x,y,z) = U_bar(x,y) * phi(z_hat) 
%                        z_hat is fractional height above bed
%       For isothermal SIA,
%        phi(z)   = (n+2)/(n+1)(1 - (1 - z_hat/H)^(n+1) )
%                    (shape function for depth variation)
%                    (see e.g. Paterson, 1994 Physics of Glaciers, p. 251)
%
%        
%      u_s(x,y)   = [(n+2/(n+1)] * U_bar(x,y)
%                    surface velocity
%   and
%      U_bar(x,y) = -(2 A)/(n+2) * (rho g)^n  ...
%                     * h(x,y)^(n+1) ( grad(S).grad(S) )^(n-1)/2 grad(S)
%
%   or, setting     C = (2 A)/(n+2) * (rho g)^n
%
%     U_bar(x,y) = - C h(x,y)^(n+1) ( grad(S).grad(S) )^(n-1)/2  grad(S) (2)
%
% Ice Flux Q(x,y)
%            Q(x,y)  = U_bar * h(x,y) 
%   can be written as
%            Q(x,y)  = alpha * U_bar * h(x,y) + (1-alpha) * U_bar * h(x,y)          (3)
%                                                             
%            Q(x,y)  = alpha * U_bar * h(x,y) - (1-alpha) * K(x,y) * grad(S)          (3)
%
%
%   Diffusion coefficient is
%            K(x,y) = - C h(x,y)^(n+2) ( grad(S).grad(S) )^(n-1)/2      (4)  
%
%   This lets us write (1) as
%
%           dh/dt + div( alpha U_bar * h(x,y) = ...
%                            div( (alpha-1) K(x,y) grad(S) ) + b_dot    (5)
%
%   which has the standard advective diffusive form -
%
%           d(rho h)/dt + div(rho u h) = div(Gamma grad(h) ) + S
%
%       For incompressible ice,    rho = rho_0
%       For advection term:     u(x,y) = alpha * U_bar(x,y)
%       For diffusion term: Gamma(x,y) = (alpha-1) * K(x,y) * rho
%       For source term:    b = rho * b_dot 
%

%  *** Complications arose using the advective term and the code has been
%  set to use only the diffusive term (s = 0)
%  Complications arise in satisfying the rule that a_P = a_E + a_W; when
%  using the upwinding scheme there is an extra term, F_e - F_w, which
%  does not go to zero in this formulation (typically sent to zero by
%  invoking mass conservation -- in this case mass conservation is our governing
%  equation!) ***

% Material Properties
%  - obeys Glen Flow Law for ice at melting point (Paterson, 1994, page 97)
%      du/dz = 2A tau^n
%      tau   = rho g (h-z) alpha           shear stress
%      A     = 6.8 10^-15 kPa^-3 s^-1      ice softness factor
%      n     = 3                           nonlinear flow law exponent
%      rho   = 900 kg m^-3                 ice density
%      g     = 9.8 m s^-2                  gravity
%  

% Boundary Conditions
%  - input flux is specified
%  - flux out the edge of the limited domain is calculated kinematically
%  (for now)
%

% Source Term
%     mass balance (net melting or accumulation) on upper surface
%
%  In the spirit of matlab, the program avoids loops over spatial indices.
%  Instead, matrices are formed for all spatially varying parameters and
%   variables.
%  Loops over time and iterates on the dependency of ice thickness at t+1
%  in the diffusivity term, where to start only ice thickness at time t is
%  available; requires accounting for underrelaxation (Patankar, CH 7)
% -----------------------------------------------------------------------