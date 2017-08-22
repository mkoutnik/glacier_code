function [ u, w, w_dzhat_dt, u_bar_P, ...
           u_bar_edges, phi, ...
           psi, int_dphi_dx, ...
           A_eff, wterm1, ...
           wterm2, wterm3, ...
           wterm4 ] = velocity_field ( ...
                                      x_grid, z_grid, zhat_grid, T, ...
                                      dS_dx, H_P, H_edges, dB_dx, b_dot, ...
                                      H_dot, B_P, S_P, d_zhat, ...
                                      A_0, x_P, x_w, x_e, ...
                                      dx_P, dx_w, dx_e, ...
                                      flux_kin, flux_kin_edges, ...
                                      W_P, W_w, W_e )
       
%--------------------------------------------------------------
%
% Michelle Koutnik (mkoutnik@ess.washington.edu)
% last updated, Sept 2007

% Velocity field from the Shallow Ice Approximation.

%  dynamic calcuation for horizontal and vertical velocity field
%  calculate u(x,z,t) and w(x,z,t) for an individual timestep.
%
%-------------------------------------------------------------



global rho_ice  g  n  Q  R  
global N_z




% integrate temperature dependence of softness parameter over depth
% ------------------------------------------------------------------

% first make a matrix, using T(x,zhat), and integrate over zhat
% zhat_grid(1,:) = 1, runs surface to bed
% (otherwise trapz and cumtrapz have problems)?

  zhat_grid_use = flipud(zhat_grid);   % now goes bed to surface in z.
  T_use         = flipud(T);
  
  
% integrate over z
% -----------------
  int_exp_QRT = cumtrapz( (exp(-Q ./ (R*T_use))) .* (1-zhat_grid_use).^n ) * d_zhat;  


% integrate over z again
% -----------------------
  int_int_exp_QRT = trapz( int_exp_QRT ) * d_zhat;  % trapz assumes dz=1.



% ----------------------
% horizontal velocity, u
% ----------------------


% dS_dx, H are the same at all depths.
% ------------------------------------
  dS_dx_grid = repmat(dS_dx, N_z, 1);
  H_grid     = repmat(H_P, N_z, 1);
  

% % u = u(x,zhat,t)
% % ---------------
    u = 2 .* A_0 .* (rho_ice * g * -dS_dx_grid).^n .* H_grid.^(n+1) .* int_exp_QRT;
   

% calculate average horizontal velocity with dynamics (SIA):
   u_bar_SIA = 2 .* A_0 .* (rho_ice * g * (-dS_dx)).^n .* H_P.^(n+1) .* int_int_exp_QRT;

   
   
   % For comparison, calculating using kinematics (ubar = Qkin / W*H) 
   u_bar_kin_compare = flux_kin ./ (W_P .* H_P);    % update so H is prescribed.
   
   u_bar_kin_edges = flux_kin_edges ./ ([W_w(1) W_e] .* H_edges);


  
% Keep dynamic calculation:
 u_bar_P = u_bar_SIA;       
 
 [ u_bar_w, u_bar_e] = get_edge_values_quadratic ( u_bar_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );
 
 u_bar_edges = [ u_bar_w(1) u_bar_e ];
   
   
 
 
 
%  figure (2)
%  plot(u_bar_P,'r')
%  hold on
%  plot(u_bar_kin_compare,'c')
 
 
 
 
 
   
% flip values so they go from surface to the bed
% -----------------------------------------------
  u = flipud(u);
  u_bar_P = flipud(u_bar_P);
  u_bar_SIA = flipud(u_bar_SIA);
  u_bar_edges = flipud(u_bar_edges);
  
   
   

% horizontal velocity shape function, phi
% ---------------------------------------
  warning('OFF', 'MATLAB:divideByZero') 
  phi = u ./ repmat(u_bar_SIA,N_z,1);
  warning('ON', 'MATLAB:divideByZero') 

  
% vertical velocity shape function, psi
% -------------------------------------
  psi = 1 - cumtrapz(phi)*d_zhat;   % psi = psi(x,zhat,t) = integral (phi)
 

  
% -----------------------
% horiztonal velocity, u -- dx/dt
% -----------------------

  u = repmat(u_bar_P,N_z,1) .* phi;     % reset here.
  
  

% --------------------
% vertical velocity, w -- dz/dt
% --------------------

  [FX] = gradient(phi);
  dphi_dx = FX ./ gradient(x_grid);   

  int_dphi_dx = cumtrapz(dphi_dx)*d_zhat;
  
  
  b_dot_grid = repmat(b_dot, N_z, 1);
  H_dot_grid = repmat(H_dot,N_z,1);
  H_grid     = repmat(H_P, N_z,1); 
 
  
  
% % -------------------  
% % OPTION 1: no melt
% % -------------------
%   mdot = 0;
%   m_dot_grid = repmat(mdot * ones(size(x_P)), N_z, 1);
% 
%   
% ----------------------
% OPTION 2: uniform melt
% ----------------------
  mdot = 0.01;   % a few cm/yr
  m_dot_grid = repmat(mdot * ones(size(x_P)), N_z, 1);
   
 
% %  % --------------------- 
% %  % OPTION 3: Patchy melt 1
% %  % ---------------------
%    melt_rate         = 0.02;
%    mdot              = zeros(size(x_P));
%    index_melt        = [2:30:length(x_P)-3];
%   for ii = 1:length(index_melt)
%        index_start = index_melt(ii); 
%        mdot(index_start:index_start+6) = melt_rate;
%    end
%    m_dot_grid        = repmat(mdot, N_z,1);
   
   
% %  % -----------------------   
% %  % OPTION 4: Patchy melt 2
% %  % -----------------------
%    melt_rate         = 0.025;
%    mdot              = zeros(size(x_P));
%    index_melt        = [2:50:length(x_P)-3];
%  
%    for ii = 1:length(index_melt)
%        index_start = index_melt(ii); 
%        mdot(index_start:index_start+20) = melt_rate;
%    end
%    m_dot_grid        = repmat(mdot, N_z,1);
   
   
   
 save mdot_use.mat m_dot_grid
   
  
 
  w = -m_dot_grid - ( (b_dot_grid - H_dot_grid - m_dot_grid) .* (psi)) ...
         + (u .* ( ((1-zhat_grid).* repmat(dB_dx,N_z,1)) ...
         +  (zhat_grid .* (dS_dx_grid)))) ...
         - ( repmat(u_bar_P,N_z,1) .* H_grid .* int_dphi_dx );


  wterm1 = m_dot_grid;
  wterm2 = ( (b_dot_grid - H_dot_grid - m_dot_grid) .* (psi));
  wterm3 = ( repmat(u_bar_P,N_z,1) .* H_grid .* int_dphi_dx );
  wterm4 = (u .* ( ((1-zhat_grid).* repmat(dB_dx,N_z,1)) + ...
           (zhat_grid .* (dS_dx_grid))));
                       

       
%  express as dzhat/dt (instead of dz/dt) for numerical integration
%  following Reeh (1988):
   w_dzhat_dt = (w ./ H_grid) - ...
               ((u./H_grid) .* ( ((1-zhat_grid).* repmat(dB_dx,N_z,1)) + ...
               (zhat_grid .* (dS_dx_grid))));


   
% so that the vector sizes are (x,z) for output
% ---------------------------------------------
  u           = u';
  w           = w';
  w_dzhat_dt  = w_dzhat_dt';
  phi         = phi';
  psi         = psi';
  int_dphi_dx = int_dphi_dx';

  

% ----------------------------------------------------
% effective isothermal softness parameter, A_eff(x,t)
% ----------------------------------------------------
  A_eff = A_0 .* (n+2) .* int_int_exp_QRT;





