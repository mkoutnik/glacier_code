

global S_in_global_Darwin

global LGM_transient


 S_in_global_Darwin   = 367;   % surf elev at the mouth of the glacier near grounding line
                  
 
% ------------
% SETUP DOMAIN
% ------------

%  set nodes in x and t defining where values
%  accumulation rate (b_dot) are defined and
%  solved for as part of the inverse problem
%  ==========================================
    [ x_nodes, t_nodes ] = load_nodes;


%  grid defining mesh for numerical solution for S(x,t)
%  ====================================================
  [ x_P, x_w, x_e,       ...
    dx_P, dx_w, dx_e,    ...
    t_P, dt_P,           ...
    z_hat, dz_hat,    ...    
    x_edges ] =   load_mesh( x_nodes, t_nodes, solver_method );


% setup size of nodal and mesh-grid values
% =========================================
    size_nodes_mesh;   % these values are GLOBAL:
                       % N_t_nodes  N_x_nodes  N_t_mesh  N_x_mesh N_z

                   
% setup grids for thermomechanical calculation
% ============================================
    x_grid       = repmat( x_P, N_z, 1 )';           
    z_hat_grid   = repmat( z_hat, N_x_mesh, 1 );  
    x_edges_grid = repmat( x_edges, N_z, 1)';
    
    
% Bed elevation B(x) -- bed does not vary in time
% ===================
      [ B_P, B_w, B_e, ...
        dB_dx_w, ...
        dB_dx_e, dB_dx_P, ...
        S_modern, ...
        value_to_zero_bed] = load_bed( x_P, x_w, x_e, dx_P, dx_w, dx_e );
      
   
% Flowband width W(x) -- width does not vary in time
% ====================
      [ W_P, W_w, W_e ] = load_width( x_P, x_w, x_e, dx_P, dx_w, dx_e );


% Enhancement factor value -- sometimes called "slip" elsewhere; sorry for that slip!
% ========================
      [ E_P, E_w, E_e, ... 
        fs_P, fs_w, fs_e ] = load_factors( x_P, x_w, x_e, dx_P, dx_w, dx_e );

    
  
% Load accumulation history b_dot(x,t)  
% ====================================
      [ b_dot_nodes, ...
        b_dot_P, ...
        b_dot_edges ] = load_b_dot( x_P, x_edges, t_P, ...
                                    x_nodes, t_nodes );
                                
      b_dot_edges_SS = b_dot_edges(1,:);

% estimate external boundary flux at left
% and right edges, as well as Q_in(t0)
% Q_ext is stacked with left values and then right values.
% =========================================
    [ Q_0_in, ...
      Q_ext_nodes, ...
      Q_ext_P ] = load_Q( t_nodes, t_P ); 

  
      

% Check if want to include tributary flux, or not.

if (add_tributary_flux == 0)
   flux_add_P = zeros(size(x_P));
end

[ flux_add_w, ...
  flux_add_e ] = get_edge_values_quadratic ...
                                  ( flux_add_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

   
                                                
  
% estimate surface elevation at x(1) and t(1)
% ===========================================
     S_0_in = S_in_global_Darwin;   % meters; defined in global_params.m

     
% estimate the surface temperature at (x_P,t)
% =========================================
    [ T_surf ] = load_T_surf; 
      

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

   T_field_0       = NaN(N_x_mesh, N_z);   % just for first time.
   T_field_0(:,1)  = T_surf(1,:);          % surface temperature boundary condition
   S_use           = repmat(S_0_in, 1, N_x_mesh);
   
  for i_x = 1:N_x_mesh
    z_use  = S_use(i_x) - B_P(i_x);
    dz_use = z_use / N_z;
    
   for i_z = 2:N_z
       T_use  = T_field_0(i_x,i_z-1);   
       K_T = 9.828 * exp( -0.0057 * T_use );  % conductivity is a function of temperature
                                              % Paterson 1994, pg. 205
       
       T_field_0(i_x,i_z) = T_use + ( (Qgeo./K_T) * dz_use ); 
   end
  end


elseif (prescribe_temperature == 1)    
      
   T_z_use = 273.15 - 3; % - 15;     
   T_field_0 = repmat(T_z_use, N_x_mesh, N_z);   % replicate everywhere 
                                  
   % Or add in something more realistic?
   

end   % if statement on temperature options
       
        

% estimate softness parameter A(T(x,t=t0))
% =========================================           
    [ A_eff_0_P, ...
      A_eff_0_w, A_eff_0_e ] = get_A_eff ( T_field_0', A_0, z_hat_grid', dz_hat, ...
                                           x_P, x_w, x_e, dx_P, dx_w, dx_e );
                                       
                                          
      
% History of surface elevation at the grounding line (GL):
% ========================================================
  S_at_GL = ones(size(t_P)) * S_0_in;   % keep constant.
   
  
  if (LGM_transient == 1)
      
    load DH_DATA/Boundary_conditions/Diamond_Hill/DH_deglaciation_scenarios.mat
    S_at_GL = S_0_in + interp1([t_P(1) -Smooth9ka.Time' 0], [Smooth9ka.HeightAboveModern(1) Smooth9ka.HeightAboveModern' Smooth9ka.HeightAboveModern(end)], t_P,'linear', 'extrap');
  
  %  S_at_GL = S_0_in + interp1([t_P(1) -Stepwise9ka.Time' 0], [Stepwise9ka.HeightAboveModern(1) Stepwise9ka.HeightAboveModern' Stepwise9ka.HeightAboveModern(end)], t_P,'linear', 'extrap');
  
    disp (' ')
    disp(' check S_at_GL in MAIN_Darwin! ')
    disp(' ')
    
  %  S_at_GL = S_0_in + interp1([t_P(1) -Smooth11ka.Time' 0], [Smooth11ka.HeightAboveModern(1) Smooth11ka.HeightAboveModern' Smooth11ka.HeightAboveModern(end)], t_P,'linear', 'extrap');
  %  S_at_GL = S_0_in + interp1([t_P(1) -Smooth15ka.Time' 0], [Smooth15ka.HeightAboveModern(1) Smooth15ka.HeightAboveModern' Smooth15ka.HeightAboveModern(end)], t_P,'linear', 'extrap');
    
  
  end



% BDM_velocity      = interp1(measures_centerline_distance, measures_flowspeed, [x_P]);
% F_check = ( (S_modern - B_P)./b_dot_P(1,:)) .* (BDM_velocity ./ x_P)
% 
% E_P = F_check / 10;
% for iii = 1:length(E_P)
%     if (E_P(iii) < 1)
%         E_P(iii) = 1;
%     end
% end
% [ E_w, E_e ] = get_edge_values_quadratic ...
%                                ( E_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );






% Initialize values
% =================
% values named "_xzt" are sized N_t_mesh, N_x_mesh, N_z  
% "_txz" would have been a better naming convention :)

T_field_xzt        = NaN( N_t_mesh, N_x_mesh, N_z );
T_field_xzt(1,:,:) = T_field_0;                       % temperature field
A_eff_edges_xt     = NaN( N_t_mesh, N_x_mesh+1 );     % effective softness parameter
A_eff_edges_xt(1,:)= [A_eff_0_w(1) A_eff_0_e];
A_eff_P_xt         = NaN( N_t_mesh, N_x_mesh );    
A_eff_P_xt(1,:)    = A_eff_0_P;
u_bar_kin_xt       = NaN( N_t_mesh, N_x_mesh+1 );     % depth-averaged horiz velocity at edges
u_bar_edges_xt     = NaN( N_t_mesh, N_x_mesh+1 );     % depth-averaged horiz velocity at edges
u_bar_P_xt         = NaN( N_t_mesh, N_x_mesh );       % depth-averaged horiz velocity at centers
phi_xzt            = NaN( N_t_mesh, N_x_mesh, N_z );  % horiz velocity shape function
psi_xzt            = NaN( N_t_mesh, N_x_mesh, N_z );  % vertical velocity shape function
int_dphi_dx_xzt    = NaN( N_t_mesh, N_x_mesh, N_z );  % integral of dphi/dx
S_P                = NaN( N_t_mesh, N_x_mesh );       % ice-surface elevation
S_w                = NaN( N_t_mesh, N_x_mesh );       % ice-surface elevation
S_e                = NaN( N_t_mesh, N_x_mesh );       % ice-surface elevation
h_P                = NaN( N_t_mesh, N_x_mesh );       % ice thickness
h_w                = NaN( N_t_mesh, N_x_mesh );       % ice thickness
h_e                = NaN( N_t_mesh, N_x_mesh );       % ice thickness
dS_dx_edges_xt     = NaN( N_t_mesh, N_x_mesh+1 );     % surface slope at edges
dS_dx_P_xt         = NaN( N_t_mesh, N_x_mesh );       % surface slope at x_P
flux_edges_dyn_xt  = NaN( N_t_mesh, N_x_mesh+1 );     % flux, values at edges
flux_edges_kin_xt  = NaN( N_t_mesh, N_x_mesh+1 );     % flux, values at edges
flux_kin_P_xt      = NaN( N_t_mesh, N_x_mesh );  
Q_out_L            = NaN( N_t_mesh, 1 );              % flux leaving end boundary
Q_out_R            = NaN( N_t_mesh, 1);
vol_change_L       = NaN( N_t_mesh, 1 );
vol_change_L(1)    = 0;
vol_change_R       = NaN( N_t_mesh, 1 );
vol_change_R(1)    = 0;
h_dot              = NaN( N_t_mesh, N_x_mesh );       % dH/dt, ice-thickness change
h_dot(1,:)         = 0;
S_dot              = NaN( N_t_mesh, N_x_mesh );       % dS/dt, surfac-elev change



 
