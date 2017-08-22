% global_variables.m


%--------------------------------------------------------------------------
%
% This script establishes and loads in global variables that are used
% throughout the code.

% (note: keep global calls to parameters that are not defined here.)

%--------------------------------------------------------------------------



% ------------------
% GENERAL PARAMETERS
% ------------------

  global itnum_max


% itnum_max   = maximum number of iterations for any loop

  itnum_max    = 50000;    

% _________________________________________________________________________ 


% ----------
% CONSTANTS
% ----------

 global s_per_year  rho_ice  g  n  Q  R
 global latent_heat_fusion
 

% s_per_year = seconds per year
% rho_ice    = density of ice [kg m^-3]
% g          = gravity [m s^-2]
% n          = flow law exponent, typically = 3
% R          = gas constant [J mol^-1 K^-1]
% latent_heat_fusion = latent heat of fusion [J/kg]
% Qgeo       = geothermal heat flux [W m^-2]

  s_per_year = 60*60*24*365.25;     
  rho_ice    = 917;                 
  g          = 9.81;               
  Q          = 60e03;             
  R          = 8.314;              

  n          = 3;   

  latent_heat_fusion = 3.34e4;  % used in initialize_temperature.m


% _________________________________________________________________________

 
% -----------------------
% FOR SURFACE CALCULATION
% -----------------------    

% global S_in_global
 global solver_method  
 global theta 
 global beta1  beta2  beta3
 global res_stop1 res_stop2 res_stop3 
 global calculate_residual

 
 
% S_in_global   = surface elevation initial condition
% solver_method = (1) implicit 1, (2) implicit 2, (3) implicit 3, (4) explicit
%                    * use solve method 2
% theta         = split between implicit/explict; 1 = fully implicit; 0 = fully explicit
% beta          = slow down convergence, for under-relaxation in implicit soln
% res_stop      = when residual is considered small enough
% reset_time    = how often to recalculate the impulse response function

% calculate_residual = flag on whether or not calculate the residual, or
%                      just iteratively update h(t+1); 1 = yes; 0 = no


 %   S_in_global   = 367;   % surf elev at the mouth of the glacier near grounding line
    
    theta         = 1;         
 
    res_stop2     = 1e-5;       % should be the only value that is used.    
    beta2         = 0.1;        % should be the only value that is used.
    
    res_stop1     = 1e-5;       % shouldn't be used                               
    res_stop3     = 6e-4;       % shouldn't be used                      
    beta1         = 0.01;       % shouldn't be used                            
    beta3         = 0.1;        % shouldn't be used 
    
    solver_method = 2;          % 1, 3 need to be debugged.

  
    calculate_residual = 0;  % 1 or 0 ; by calculating the residual it takes about
                             % six+ times longer! the same answer is
                             % achieved with either calculation. so, keep = 0.

   

% ---------------------------------
% For Impulse Response Calculation
% ---------------------------------        
    
  global reset_time
  global use_vialov_ext
  global paterson_c_over_a
  
  global dt_impulse
  
  global ss_return_value

  global perturbation_vol
  
  
% reset_time  = how often to recalculate impulse response function    
% use_vialov_ext        = which shape to extend limited model with - flag. 
%                         (1) vialov; (0) paterson
% paterson_c_over_a     = ratio of accumulation rate to ablation rate 
%                         for paterson extension of limited domain
% dt_impulse            = timestep for calculation of full domain response
%                         to impulse perturbation
% perturbation_vol      = volume addedfor impulse response calculation


    reset_time          = 2000;       % years

    use_vialov_ext      = 0;
    
    paterson_c_over_a   = 1;           % 1;

    dt_impulse          = 10;          % timestep in years
                                       % used in calc_impulse_response.m
    
    
    ss_return_value     = 1e-6;
      
    perturbation_vol    = 10;         % always add the same volume
    
                             
                             
% % _________________________________________________________________________  
%   
% % ------------------------
% % FOR THERMAL CALCULATION
% % ------------------------
% 
% % in thermomechanical_calc.m:
%   global max_T_change  
% 
% % in make_domain.m:
%   global zscale xscale
%   global z_thin_vols
%   global nc_start nr_start
% 
% % in initialize_temperature.m:
%   global K_rock c_rock rho_rock
%   global k_ini_ice
% 
% % in thermal_main.m:
%   global thermal_flag
%   global K_ice c_ice
% 
% 
% % ** additional flags and values set in curvilinear_cs.m **
% 
% 
% % max_T_change = temperature difference between iterations, this loop makes
% %                the code "thermomechanical"
% % zscale       = length scale used to "stretch" vertical unit vectors for grid spacing
% % xscale       = length scale used to "stretch" horizontal unit vectors for grid spacing
% % z_thin_vols  = "=1" to include a thin volume at bed for melting rate calculation 
% %                "=0" for no thin volumes at bed,
% % nc_start     = nc - 1 is the number of volumes in the x direction 
% %               ( nc is  the number of volume faces in the x direction )
% % nr_start     =  nr - 1 is the number of volumes in the z direction 
% %               ( nr is  the number of volume faces in the z direction )
% 
% % if thermal_flag = 1, temperature dependent thermal properties in ice
% %                 = 0, constant thermal properties in ice       
% 
% 
%   zscale       = 1000;        
%   xscale       = 40000;              
%   z_thin_vols  = 0;         
%   nc_start     = 80;                 
%   nr_start     = 60;                
% 
%   max_T_change = 0.02;  % degrees
% 
%   K_rock       = 2.0;
%   c_rock       = 1000; 
%   rho_rock     = 2700;
%   k_ini_ice    = 6.6e7;
% 
%   thermal_flag = 1
% 
%   K_ice        = 6.6e7;
%   c_ice        = 2093;
% 
% 
% % _________________________________________________________________________
% 
% 
% % ----------------------  
% % FOR PARTICLE TRACKING
% % ----------------------
% 
% global dx_back_off  N_paths  N_paths_track
% global T_span_interval
% global N_paths_make_layers
% 
% global h_div  h_flank
% 
% 
% % dx_back_off = to prevent any paths exiting the domain during iterations, downstream 
% %               path should reach layer at dx_back_off fraction of x_node domain
% % N_paths     = number of particle paths to track
% % T_span_interval = used in particle_main.m for ode113 solver for particle
% %                   tracking
% % h_div       = 'h' for divide in Dansgaard-Johnsen model. See Nereson and Waddington 2002
% % h_flank     = 'h' for flank in Dansgaard-Johnsen model. See Nereson and Waddington 2002
% % N_paths     = the number of particle paths tracked; layer-depth data is
% %               interpolated to the positions of the modelled layer
% 
% 
%   dx_back_off         = 0.05;
%   
%   N_paths              = 500; 
%   
%   N_paths_track        = 500;
%   
%   N_paths_make_layers = 500;
%   
%   T_span_interval     = 500;
%     
%   h_div               = 0.6;
%   h_flank             = 0.2;
%  
% % _________________________________________________________________________
% 
% 
% % ---------------------  
% % FOR INVERSE PROBLEM
% % ---------------------
% 
% 
% global N_data
% global N_u N_b N_S
% global N_layers
% global M_param
% global N_h_data
% global N_core
% % declare global here, but assign in MAIN.
% 
% global stop_test
% global N_segments
% global solve_with_svd
% global num_perturb_factor
% 
% 
% % stop_test          = when to stop iterating on nu value
% % N_segments         = number of intergration intervals for numerical integration
% %                    of du_i/dp_j in dxs_dp.m
% % solve_with_pinv    = use pinv or SVD to solve for changes in parameters
% %                      1 = use pinv; 0 = use SVD
% % num_perturb_factor = percent by which to perturb values for calculation
% %                      of numerical jacobians
% 
% 
%  stop_test          = 1000;
%  N_segments         = 500;    
%  solve_with_svd     = 1; 
%  
 
 
