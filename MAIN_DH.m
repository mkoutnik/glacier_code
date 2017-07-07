
% Transient glacier surface evolution model
% -----------------------------------------
% close all; 
% Last updated April 2017
% Michelle Koutnik

% Using the finite (control) volume method (FVM)
% Shallow Ice Approximation with optional parameterizations for:
%       - sliding (simple partitioning; Cuffey and Paterson pg. XX)
%       - lateral drag (Adhikari and Marshall 2012)
%       - longitudinal stress gradients (e.g., Kamb and Echelmeyer 1986;
%       Adhikari and Marshall 2011)

% USER DEFINED:
% Input data from specified folder (e.g., /BDM_DATA)
% Files with boundary conditions from specified folder (e.g., /TEST_BDM)

% In MAIN.m: run flags; global variables; heat flux; softness parameter;
% temperature as a function of depth; ice-surface elevation initial
% condition at the grounding line

% Primary boundary condition to change is surface elevation at the
% grounding line ("S_at_GL"); could also change flux boundary conditions at
% upstream and downstream ends of the domain ("Q_ext"); 
% could also change surface mass balance ("b_dot")

% Options available to use simple minimization search to better constrain poorly
% known parameter values / boundary conditions


% -------------------------------------------------------------------------

addpath('Factorize')             % use fast routines for matrix inversion.
addpath('regularization_tools')  % tools for numerical regularization.
                                 % http://www2.imm.dtu.dk/~pch/Regutools/
disp(' ')
disp(' ')

global DIRECTORY_data DIRECTORY_setup

global linear_temperature prescribe_temperature

global min_search_E min_search_fs min_search_bed

global deformation_only deformation_plus_sliding sliding_only ...
%       deformation_sliding_lateraldrag deformation_sliding_longstress ...
%       deformation_sliding_lateraldrag_longstress

global figures_on

global lower_resolution

% Simple way to try and find values of boundary conditions / parameters 
% that may be poorly constrained:
% "E": enhancement factor on shear deformation
% "bed": bed topography
% "width": flowband width
% "fs": sliding factor in simple partitioning between deformation/sliding

% Refine estimates of these parameters - one at a time! - to better match
% ice-surface elevation and/or ice-surface velocity data
% Need to set steady_state_only = 1 when running min search
min_search_E     = 0;   
min_search_fs    = 0;   
min_search_bed   = 0;

lower_resolution = 0;   % Runs faster. Use spatial step of multiple km.

figures_on       = 0;   % to set this =1 must have run once with all min_search=0,
                        % then can set=1 with one of the min_search=1;
                        % plots compare difference between two runs

                        
steady_state_only = 1;  % Flag in loop below so that only run one 
                        % calculation to get steady state values.
                        % *Also if want to run minimization search*
                        
                        
% Only one of these should be set = 1 at a time. Compare output from
% different assumptions about ice flow
deformation_only                           = 1;
deformation_plus_sliding                   = 0; 
sliding_only                               = 0; 
deformation_sliding_lateraldrag            = 0; % not included yet!
deformation_sliding_longstress             = 0; % not included yet!
deformation_sliding_lateraldrag_longstress = 0; % not included yet!


% initial value, or initial *guess* of temperature: must set one of these = 1!      
linear_temperature       = 0;   % 1 = temperature varies linearly with depth
                                % 0 = don't use this formulation.                                
prescribe_temperature    = 1;   % 1 = prescribe how temperature varies with depth
                                % 0 = don't use this formulation.



% Directory pointers
% ------------------
  DIRECTORY_data = 'DH_DATA';
  addpath( DIRECTORY_data )
  DIRECTORY_setup = 'TEST_DH'; 
  addpath( DIRECTORY_setup )

                    
              
% Darwin-Hatherton data
% ================
load DH_accum_width_velocity.mat
load DH_tributaries.mat



% load all other globals:       
% -----------------------
  global_variables;


                                
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

      
% estimate surface elevation at x(1) and t(1)
% ===========================================
     S_0_in = S_in_global;   % meters; defined in global_params.m

     
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
      
   T_z_use = 273.15 - 25; % - 15;     
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
%  S_at_GL = chronology_data (t_P);      % use data.
   



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



 

% -------------------------------------------------------------------------
if ( (min_search_E == 1) || (min_search_bed == 1) || (min_search_fs == 1) )
    
  time = 1;  
    
  if (min_search_E == 1)
    run_min_search_E;
  end   % end on min_search

  if (min_search_fs == 1)
    run_min_search_fs;
  end   % end on min_search

  if (min_search_bed == 1)
    run_min_search_bed;
  end   % end on min_search

end
% -------------------------------------------------------------------------




disp(' ' )  
disp('    Ice-surface calculation...')
disp(' ' )

if (steady_state_only == 1)

    time = 1;
    S_at_GL(time) = S_0_in;  % Reset since t_P(1) is usually 20ka, not present
    disp('Steady state only')
    
    calculate_surface;
    
    
else
    
 for time = 1:N_t_mesh
     
   if (time == 1)
       disp('Start in steady state')
   end
   if ( (time == 1) && (min_search_bed == 1) || (min_search_E == 1) || (min_search_fs == 1))
       disp('Start in steady state with minimized parameters')
   end
       disp(['Time in calculation = ',int2str(t_P(time))]);
       
       calculate_surface;         % script

 end   % for loop on time

end  % if statement on steady_state


                 
save save_output.mat





figure(10)   % SURFACE and BED
plot(x_P/1000, B_P,'r', 'linewidth', 2)
hold on
plot(x_P/1000, S_modern, 'b', 'linewidth', 2)
plot(x_P/1000, S_P(1,:), 'linewidth', 2)
legend('Bed', 'Measured surface', 'Calculated surface', 'location', 'northwest')
title('Surface and Bed topography')
xlabel('Distance along flowband (km)')
ylabel('Elevation (m)')

figure(20)   % ACCUMULATION
hold on
plot(x_P/1000, b_dot_P', 'b')
title('Accumulation Rate')
xlabel('Distance along flowband (km)')
ylabel('Accumulation rate (m/yr)')

figure(30)   % SURFACE VELOCITY 
plot(x_edges/1000, average_vel_estimate*(5/4), 'k')
hold on
plot(Darwin_measures_centerline_distance/1000, Darwin_measures_flowspeed,'m', 'linewidth', 1)
title('Surface velocity')
xlabel('Distance along flowband (km)')
ylabel('Velocity (m/yr)')
legend('Calculated velocity',  'MEaSUREs velocity', 'location', 'best')

% figure(40)
% plot(t_P/1000, S_at_GL-S_at_GL(end), 'b')
% title('Change in elevation at grounding line since present day')
% ylabel('Elevation (m)')
% xlabel('Time (kyr)')
% ylim([-5 800])
%  


