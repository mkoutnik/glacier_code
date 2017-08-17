


global S_in_global_Hatherton


 S_in_global_Hatherton   = 793; 
                  
                                
% ------------
% SETUP DOMAIN
% ------------

%  set nodes in x and t defining where values
%  accumulation rate (b_dot) are defined and
%  solved for as part of the inverse problem
%  ==========================================
    [ x_nodes2, t_nodes2 ] = load_nodes2;


%  grid defining mesh for numerical solution for S(x,t)
%  ====================================================
  [ x_P2, x_w2, x_e2,       ...
    dx_P2, dx_w2, dx_e2,    ...
    t_P2, dt_P2,           ...
    z_hat2, dz_hat2,    ...    
    x_edges2 ] =   load_mesh2( x_nodes2, t_nodes2, solver_method );


% setup size of nodal and mesh-grid values
% =========================================
    size_nodes_mesh2;   % these values are GLOBAL:
                       % N_t_nodes  N_x_nodes  N_t_mesh  N_x_mesh N_z

                   
% setup grids for thermomechanical calculation
% ============================================
    x_grid2       = repmat( x_P2, N_z2, 1 )';           
    z_hat_grid2   = repmat( z_hat2, N_x_mesh2, 1 );  
    x_edges_grid2 = repmat( x_edges2, N_z2, 1)';
    
    
% Bed elevation B(x) -- bed does not vary in time
% ===================
      [ B_P2, B_w2, B_e2, ...
        dB_dx_w2, ...
        dB_dx_e2, dB_dx_P2, ...
        S_modern2, ...
        value_to_zero_bed2] = load_bed2( x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );
      
   
% Flowband width W(x) -- width does not vary in time
% ====================
      [ W_P2, W_w2, W_e2 ] = load_width2( x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );


% Enhancement factor value -- sometimes called "slip" elsewhere; sorry for that slip!
% ========================
      [ E_P2, E_w2, E_e2, ... 
        fs_P2, fs_w2, fs_e2 ] = load_factors2( x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );

    
  
% Load accumulation history b_dot(x,t)  
% ====================================
      [ b_dot_nodes2, ...
        b_dot_P2, ...
        b_dot_edges2 ] = load_b_dot2( x_P2, x_edges2, t_P2, ...
                                    x_nodes2, t_nodes2 );
                                
      b_dot_edges_SS2 = b_dot_edges2(1,:);

% estimate external boundary flux at left
% and right edges, as well as Q_in(t0)
% Q_ext is stacked with left values and then right values.
% =========================================
    [ Q_0_in2, ...
      Q_ext_nodes2, ...
      Q_ext_P2 ] = load_Q2( t_nodes2, t_P2 ); 

  
   
  
  
if (add_tributary_flux == 0)
   flux_add_P2 = zeros(size(x_P2));
end

[ flux_add_w2, ...
  flux_add_e2 ] = get_edge_values_quadratic ...
                                  ( flux_add_P2, x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );

   
                                               
     
% estimate the surface temperature at (x_P,t)
% =========================================
    [ T_surf2 ] = load_T_surf2; 
      

% Estimate the heat flux
% =======================
    Qgeo = 0.065;  % W m-2
    
    
% estimate global parameter A_0
% =============================
    A_0 = 8.0e-04;  % kPa^(-3) s^(-1); Paterson (1994, pg. 97)
    A_0 = A_0 * (s_per_year / 1.0e09);  % convert to Pa-3 yr-1      
   
   

    
% temperature field T(x,z)
% ========================    
if (linear_temperature == 1)   % temperature varies linearly from surface to bed.

%    T_field_0       = NaN(N_x_mesh, N_z);   % just for first time.
%    T_field_0(:,1)  = T_surf(1,:);          % surface temperature boundary condition
%    S_use           = repmat(S_0_in, 1, N_x_mesh);
%    
%   for i_x = 1:N_x_mesh
%     z_use  = S_use(i_x) - B_P(i_x);
%     dz_use = z_use / N_z;
%     
%    for i_z = 2:N_z
%        T_use  = T_field_0(i_x,i_z-1);   
%        K_T = 9.828 * exp( -0.0057 * T_use );  % conductivity is a function of temperature
%                                               % Paterson 1994, pg. 205
%        
%        T_field_0(i_x,i_z) = T_use + ( (Qgeo./K_T) * dz_use ); 
%    end
%   end


elseif (prescribe_temperature == 1)    
      
   T_z_use2 = 273.15 - 5; % - 15;     
   T_field_02 = repmat(T_z_use2, N_x_mesh2, N_z2);   % replicate everywhere 
                                  
   % Or add in something more realistic?
   

end   % if statement on temperature options
       
        

% estimate softness parameter A(T(x,t=t0))
% =========================================           
    [ A_eff_0_P2, ...
      A_eff_0_w2, A_eff_0_e2 ] = get_A_eff ( T_field_02', A_0, z_hat_grid2', dz_hat2, ...
                                           x_P2, x_w2, x_e2, dx_P2, dx_w2, dx_e2 );
                                       
% estimate surface elevation at x(1) and t(1)
% ===========================================
     S_0_in2 = S_in_global_Hatherton;   % meters; defined in global_params.m
                                          
      
% History of surface elevation at the grounding line (GL):
% ========================================================
  S_at_GL2 = ones(size(t_P2)) * S_0_in2;   % keep constant / first estimate
   




% Initialize values
% =================
% values named "_xzt" are sized N_t_mesh, N_x_mesh, N_z  
% "_txz" would have been a better naming convention :)

T_field_xzt2        = NaN( N_t_mesh2, N_x_mesh2, N_z2 );
T_field_xzt2(1,:,:) = T_field_02;                       % temperature field
A_eff_edges_xt2    = NaN( N_t_mesh2, N_x_mesh2+1 );     % effective softness parameter
A_eff_edges_xt2(1,:)= [A_eff_0_w2(1) A_eff_0_e2];
A_eff_P_xt2         = NaN( N_t_mesh2, N_x_mesh2 );    
A_eff_P_xt2(1,:)    = A_eff_0_P2;
u_bar_kin_xt2       = NaN( N_t_mesh2, N_x_mesh2+1 );     % depth-averaged horiz velocity at edges
u_bar_edges_xt2     = NaN( N_t_mesh2, N_x_mesh2+1 );     % depth-averaged horiz velocity at edges
u_bar_P_xt2         = NaN( N_t_mesh2, N_x_mesh2 );       % depth-averaged horiz velocity at centers
phi_xzt2            = NaN( N_t_mesh2, N_x_mesh2, N_z2 );  % horiz velocity shape function
psi_xzt2            = NaN( N_t_mesh2, N_x_mesh2, N_z2 );  % vertical velocity shape function
int_dphi_dx_xzt2    = NaN( N_t_mesh2, N_x_mesh2, N_z2 );  % integral of dphi/dx
S_P2                = NaN( N_t_mesh2, N_x_mesh2 );       % ice-surface elevation
S_w2                = NaN( N_t_mesh2, N_x_mesh2 );       % ice-surface elevation
S_e2                = NaN( N_t_mesh2, N_x_mesh2 );       % ice-surface elevation
h_P2                = NaN( N_t_mesh2, N_x_mesh2 );       % ice thickness
h_w2                = NaN( N_t_mesh2, N_x_mesh2 );       % ice thickness
h_e2                = NaN( N_t_mesh2, N_x_mesh2 );       % ice thickness
dS_dx_edges_xt2     = NaN( N_t_mesh2, N_x_mesh2+1 );     % surface slope at edges
dS_dx_P_xt2         = NaN( N_t_mesh2, N_x_mesh2 );       % surface slope at x_P
flux_edges_dyn_xt2  = NaN( N_t_mesh2, N_x_mesh2+1 );     % flux, values at edges
flux_edges_kin_xt2  = NaN( N_t_mesh2, N_x_mesh2+1 );     % flux, values at edges
flux_kin_P_xt2      = NaN( N_t_mesh2, N_x_mesh2 );  
Q_out_L2            = NaN( N_t_mesh2, 1 );              % flux leaving end boundary
Q_out_R2            = NaN( N_t_mesh2, 1);
vol_change_L2       = NaN( N_t_mesh2, 1 );
vol_change_L2(1)    = 0;
vol_change_R2       = NaN( N_t_mesh2, 1 );
vol_change_R2(1)    = 0;
h_dot2              = NaN( N_t_mesh2, N_x_mesh2 );       % dH/dt, ice-thickness change
h_dot2(1,:)         = 0;
S_dot2              = NaN( N_t_mesh2, N_x_mesh2 );       % dS/dt, surfac-elev change



 


