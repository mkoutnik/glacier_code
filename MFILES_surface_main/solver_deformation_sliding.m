function [ S_P_t_guess, ...
           S_P_t, ...
           flux_edges_dyn_t ] = solver_deformation_sliding( flow_constant, flow_constant_t, ...
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
                                   deformation_only, deformation_plus_sliding, sliding_only)
                               
                               
% -------------------------------------------------------------------------            
                              
global n
global theta
global beta2
global rho_ice g 
   


% set up interface "conductivity"
% ===============================

% Can't have negative coefficients so take (dS_dx_e.^2).^( (n-1)/2 )

    K_w   = ( (( E_w .* flow_constant(1:end-1) .* W_w .* (h_w.^(n+2))) ...
                + (fs_w .* ((rho_ice*g)^n .* h_w.^(n) .* W_w)) ) ...
            .* ( (dS_dx_w.^2).^( (n-1)/2 ) ) );
        
    K_e   = ( (( E_e .* flow_constant(1:end-1) .* W_e .* (h_e.^(n+2))) ...
                + (fs_e .* ((rho_ice*g)^n .* h_e.^(n) .* W_e)) ) ...
            .* ( (dS_dx_e.^2).^( (n-1)/2 ) ) );     
   
    K_w_t = ( (( E_w .* flow_constant(1:end-1) .* W_w .* (h_w_t.^(n+2))) ...
                + (fs_w .* ((rho_ice*g)^n .* h_w_t.^(n) .* W_w)) ) ...
            .* ( (dS_dx_w_t.^2).^( (n-1)/2 ) ) );    

    K_e_t = ( (( E_e .* flow_constant(1:end-1) .* W_e .* (h_e_t.^(n+2))) ...
                + (fs_e .* ((rho_ice*g)^n .* h_e_t.^(n) .* W_e)) ) ...
            .* ( (dS_dx_e_t.^2).^( (n-1)/2 ) ) );        

      
% set up interface "diffusion coefficient"
% ========================================
     D_w   = (K_w) .* 1./dx_w;   % units are m^2/yr
     D_e   = (K_e) .* 1./dx_e;

     D_w_t = (K_w_t) .* 1./dx_w; 
     D_e_t = (K_e_t) .* 1./dx_e;
 
     
     
     
% use solution algorithm from Patankar (1980; pg. 54-)
% form coefficients
% =================   
     a_W       = theta * (1./W_P) .* D_w_t;    % upstream FV centerpoint
     a_E       = theta * (1./W_P) .* D_e_t;    % downstream FV centerpoint

     a_W(1)    = 0;    % boundary condition -- thickness
     a_E(end)  = 0;    % boundary condition -- flux
    
     dx_bc_1   = 0.5 * dx_P(1);               % for boundary condition, 
     dx_bc_end = 0.5 * dx_P(end);             % NOT using half volumes. ?
    
     
     a_P0      = (dx_P./dt_P);                % coefficient for time-1
     
     a_P0(1)   = dx_bc_1 / dt_P(1);   % boundary condition over half volume?
     
     a_P       = a_P0 + a_E + a_W;            % FV centerpoint
     
     a_P(1)    = a_P(1) + 2*(theta * (1./W_P(1)) .* D_w_t(1));
     
% Now, calculate the flux dynamically.
% ====================================
     Q_w      = - sign(dS_dx_w) .* K_w .* abs(dS_dx_w);   % flux at western edges at time=time-1
     Q_e      = - sign(dS_dx_e) .* K_e .* abs(dS_dx_e);   % flux at eastern edges at time=time-1

     Q_w_t    = - sign(dS_dx_w_t) .* K_w_t .* abs(dS_dx_w_t);   % flux at western edges at time=time
     Q_e_t    = - sign(dS_dx_e_t) .* K_e_t .* abs(dS_dx_e_t);   % flux at eastern edges at time=time


% check that this is the same as dynamic calculation: 
  A_eff_edges_t = flow_constant_t * ((n+2)/2) * (1/((rho_ice*g)^n));
  Q_e_dyn = calc_flux_dyn( x_e, h_e_t, dS_dx_e_t, E_e, fs_e, W_e, ...
                           A_eff_edges_t(2:end), deformation_only, deformation_plus_sliding, sliding_only );  
  
%   plot(S_e_t)
%   hold on
%   
%   plot(x_e, Q_e_dyn,'r')
%   hold on
%   plot(x_e, Q_e_t,'c--')
%                        
% _________________________________________________________________________

% Dynamic:
% -----------
  flux_edges_dyn_t = [Q_w_t(1) Q_e_t];

% _________________________________________________________________________                          
 
                                    

% form right-hand-side vector of known values
% ===========================================
     rhs      = (S_P .* a_P0) - ...
                ( (1./W_P) .* (1-theta) .* (Q_e - Q_w) ) + ...
                (( (theta * b_dot_P(2,:)) + ...
                ((1-theta) * b_dot_P(1,:)) ) .* dx_P);
            
% Prescribe surface elevation at the grounding line:
   %  rhs(1)    = rhs(1) + (2*(S_at_GL-B_P(1)) * theta * (1./W_P(1)) .* D_w_t(1));
   
     rhs(1)    = rhs(1) + (2*(S_at_GL(1)) * theta * (1./W_P(1)) .* D_w_t(1));
    
    
% Prescribe flux at head = SS value (nothing added):
%     flux_in   = Q_out_R_SS;  
%     flux_in_t = Q_out_R_SS;  
     flux_in   = Q_out_R_SS - Q_external_R(1);  
     flux_in_t = Q_out_R_SS - Q_external_R(2); 
     
     rhs(end) = (S_P(end) .* a_P0(end)) - ...   % Or, h_P(end)?
                ( (1/W_P(end))  * theta * flux_in_t ) - ...
                ( (1/W_P(end)) * (1-theta) * (flux_in - Q_w(end)) )  + ...
                (( (theta * b_dot_P(2,end)) + ...
                 ((1-theta) * b_dot_P(1,end)) ) .* dx_bc_end);
             
    

  S_P_t_guess = S_P_t;   % guess h(t) = h(t-1); store this value.


     
% solve
% =====
    M = diag(a_P) + diag(-a_E(1:N_x-1),1) + diag(-a_W(2:N_x),-1);
    % tridiagonl matrix.
 
%     soln_pinv = M\rhs';   % slower!   
    
    F_temp    = factorize(M);   % use files from matlab central.
    soln_pinv = F_temp\rhs';    % must addpath('Factorize') to access.
    
    
    S_P_t_soln = soln_pinv';   % new estimate of S_P_t
    

    
% slow down convergence in order for stability of iterative solution!
% must invoke underrelaxation
% ============================
  S_P_t = ((1-beta2)*S_P_t_guess)+(beta2*S_P_t_soln);
  
                    

   
   
                             