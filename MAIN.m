% MAIN -- calls MAIN_Darwin and MAIN_Hatherton

% First is to generate surface for Darwin, then calculate surface for
% Hatherton based on matching to Darwin at 17 km along flowline, then
% recalculate for Darwin with flux from Hatherton added 



% Transient glacier surface evolution model
% -----------------------------------------

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

global min_search_E min_search_fs min_search_bed min_search_E_and_fs

global deformation_only deformation_plus_sliding sliding_only ...
%       deformation_sliding_lateraldrag deformation_sliding_longstress ...
%       deformation_sliding_lateraldrag_longstress

global deformation_only2 deformation_plus_sliding2 sliding_only2 ...
%       deformation_sliding_lateraldrag2 deformation_sliding_longstress2 ...
%       deformation_sliding_lateraldrag_longstress2

global add_tributary_flux

global figures_on

global lower_resolution

global steady_state_only

global LGM_steady_state

global LGM_transient


% - - - - - - - - -  - - - - - - - - - - - - -  

% -------------------------------------------------------------------------
% STEP 1: Run minimization scheme to find best E and fs that give match to
% surface elevation and surface velocity at Darwin and at Hatherton

steady_state_only = 0;  % =1 to run.

% Need to set one of the min_search algorithms = 1 at a time.
% THEN, need to set corresponding flags for deformation and/or sliding.

% Simple way to try and find values of boundary conditions / parameters 
% that may be poorly constrained:
% "E": enhancement factor on shear deformation
% "bed": bed topography
% "width": flowband width
% "fs": sliding factor in simple partitioning between deformation/sliding
% Refine estimates of these parameters - one at a time! - to better match
% ice-surface elevation and/or ice-surface velocity data
min_search_E        = 0;   
min_search_fs       = 0; 
min_search_E_and_fs = 0;
min_search_bed      = 0;  % not used!

lower_resolution = 0;   % Runs faster. Use spatial step of multiple km.
                        % Likely not as accurate, especially for gradients
                        % at boundaries (e.g., dS/dx)

                        
% - - - - - - - - - - - - - - - - - - - - - - -                         


% - - - - - - - - - - - - - - - - - - - - - - -                                                 
% Only one of these should be set = 1 at a time. 
% For Darwin:
deformation_only                           = 0;
sliding_only                               = 0; 
deformation_plus_sliding                   = 1; 
deformation_sliding_lateraldrag            = 0; % not included yet!
deformation_sliding_longstress             = 0; % not included yet!
deformation_sliding_lateraldrag_longstress = 0; % not included yet!

% For Hatherton (always as XX2)
deformation_only2                           = 0;
sliding_only2                               = 0; 
deformation_plus_sliding2                   = 1; 
deformation_sliding_lateraldrag2            = 0; % not included yet!
deformation_sliding_longstress2             = 0; % not included yet!
deformation_sliding_lateraldrag_longstress2 = 0; % not included yet!
% - - - - - - - - - - - - - - - - - - - - - - -                         


% Go into run_min_search_E, run_min_search_fs, or run_min_search_E_and_fs
% and run separately for Darwin (1) and Hatherton (2) need to set at top of iteration
% loop. I should make this better, sorry -- if I get time to do it, I will!

% After each minimization run I saved one file for each glacier:
% minE_D.mat: minimization for E for Darwin
% minE_H.mat: minimization for E for Hatherton
% minfs_D.mat: minimization for fs for Darwin
% minfs_H.mat: minimization for fs for Hatherton
% minE_and_fs_D.mat: minimization for E and fs for Darwin
% minE_and_fs_H.mat: minimization for E and fs for Hatherton

% After have all of these files loaded then run snipet of code in
% set_boundary conditions to generate three files that have values for
% Darwin and Hatherton, and that will be loaded in step2.

% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% STEP 2: Run these cases with prescribed boundary conditions to get steady
% state surface to compare to modern data

% Set steady_state_only = 1 (above)

% Go back up and change flags for deformation / sliding to be correct pair!

% For cases 1-3 need to uncomment snipet in set_boundary_conditions
% corresponding to E and fs patterns that want to use (those found from
% minization)


% As of now, there are 6 cases that I tried:
% Case 1. deformation_only = 1; deformation_only2 = 1 and min_E_values in set_boundary_conditions.m
%    Run saved as min1.mat
% Case 2. sliding_only = 1; sliding_only2 = 1 and min_fs_values in set_boundary_conditions.m
%    Run saved as min2.mat
% Case 3. deformation_plus_sliding = 1; deformation_plus_sliding2 = 1 and minE_and_fs_values.mat
%    Run saved as min3.mat
% Case 4. deformation_only = 1; deformation_only2 =1 
%    Run saved as min4.mat
% Case 5. sliding_only = 1; sliding_only2 = 1
%    Run saved as min5.mat
% Case 6. sliding_only = 1 (Darwin); deformation_only2 = 1 (Hat)
%     Run saved as min6.mat
% Case 7. deformation_only = 1; sliding_only2 = 1
%     Run saved as min7.mat   -- not always carried around...

% Could try more of these, or different combinations of deformation or
% sliding for Darwin or Hatherton. Didn't get there yet, but eventually
% could.

% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% STEP 3: Run these cases with prescribed boundary conditions to get steady

% Set run_min... = 0 and keep steady_state = 1; set LGM_steady_state=1                        
  

LGM_steady_state = 0;   % Run for steady state to best match LGM limits on Hatherton
                        % Set only this steady state flag if want to do this.
                        % Probably want lower_resolution = 0 for this one.
     
% Run this for each case separately and set boundary conditions appropriately, 
% same as above. Saved each as LGM_steady_state_min1...6                        
% -------------------------------------------------------------------------                        
                        
% -------------------------------------------------------------------------
% STEP 4: Run these cases with prescribed boundary conditions to get
% transient
  
% Set run_min... = 0 and keep steady_state = 0; set LGM_steady_state = 0
% set LGM_transient = 1

% Pick set of boundary conditions to use for the run (following cases
% above, or another)

% Need to prescribe elevation at mouth of Darwin, this is set in
% MAIN_Darwin at the bottom as "S_at_GL" -- pick between different smooth
% and step histories to compare to data
                        
LGM_transient    = 1;   % Sets prescribed S_at_GL for Darwin



% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Make sure all other boundary conditions are set the way you want them,
% basically everything in TEST_DH directory. 

% In particular:
% load_b_dot, load_b_dot2
% load_width, load_width2

% Figures can auto generate with these output files in plot_figures.m

% -------------------------------------------------------------------------



add_tributary_flux       = 0;   % add flux from tributaries at locations 
                                % along length of flowline; keep=0 for now.

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
  
  

  
% -------------------------------------------------------------------------
  MAIN_Darwin
  disp('Loading for Darwin')
% -------------------------------------------------------------------------


  xpos_Hatherton     = 35000;  % Where modern surface elevations same from data; grid starts at 10 km
  index_xpos_Hatherton_P = find (x_P <= xpos_Hatherton, 1, 'last');
  index_xpos_Hatherton_edges = find (x_w <= xpos_Hatherton, 1, 'last');
  

% -------------------------------------------------------------------------
  MAIN_Hatherton
  disp('Loading for Hatherton')
% -------------------------------------------------------------------------




disp(' ' )  
disp('    Ice-surface calculation...')
disp(' ' )



if (steady_state_only == 1)

     
% -------------------------------------------------------------------------
if ( (min_search_E == 1) || (min_search_bed == 1) || (min_search_fs == 1) || (min_search_E_and_fs == 1) )
    
  time  = 1;  
  time2 = 1; 
  
  if (min_search_E_and_fs == 1)
    run_min_search_E_and_fs;    % run_min_search_E_2;
  end  
  
  if (min_search_E == 1)
    run_min_search_E;    % run_min_search_E_2;
  end   % end on min_search

  if (min_search_fs == 1)
    run_min_search_fs;
  end   % end on min_search

  if (min_search_bed == 1)
    run_min_search_bed;
  end   % end on min_search

% -------------------------------------------------------------------------
    
    
elseif (LGM_steady_state == 1)   % find LGM thickness at Darwin GL that 
                                 % best matches Hatherton limits

  time  = 1;  
  time2 = 1; 

% -------------------------------------------------------------------------  
%   Set different boundary conditions from minimization runs -- 
%     reset values for E, fs, and deformation/sliding flag!
  set_boundary_conditions_LGM;   % for now, manually pick which one!
% -------------------------------------------------------------------------  
  
  
% -------------------------------------------------------------------------  
  run_min_search_GL_at_LGM;
% -------------------------------------------------------------------------
   
   


   
else   % run for steady state and/or transient with no minimization.    
    
    time  = 1;
    time2 = 1;
    
    S_at_GL(time) = S_0_in;  % Reset since t_P(1) is usually 20ka, not present
    
    
 % -------------------------------------------------------------------------         
    set_boundary_conditions_LGM;   % use this file for now. 
% -------------------------------------------------------------------------  
    
    
    disp('Steady state only.')
    disp(' ')
    
% -------------------------------------------------------------------------    
    calculate_surface_Darwin;
% -------------------------------------------------------------------------

      disp('Calculating for Darwin')

      
% RESET Hatherton from Darwin calculation
disp('  Resetting surface elevation for Hatherton based on Darwin values ... ')
%S_in_global_Hatherton = h_w(time2, index_xpos_Hatherton_edges) + B_w(index_xpos_Hatherton_edges) - value_to_zero_bed + value_to_zero_bed2;
S_in_global_Hatherton = h_w(time2, index_xpos_Hatherton_edges) + B_w(index_xpos_Hatherton_edges);
S_0_in2               = S_in_global_Hatherton;
S_at_GL2              = ones(size(t_P2)) * S_0_in2;  
      
     
% -------------------------------------------------------------------------
    calculate_surface_Hatherton;
      disp('Calculating for Hatherton')
% -------------------------------------------------------------------------
    
    
% Then, setup flux added at 17 km from Hatherton...
flux_add_P(index_xpos_Hatherton_P) = flux_edges_dyn_xt2(time, 1) / W_P2(1);
[ flux_add_w, ...
  flux_add_e ] = get_edge_values_quadratic ...
                                  ( flux_add_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );                              
                              
% -------------------------------------------------------------------------                              
   calculate_surface_Darwin;
      disp('Recalculating for Darwin')
% -------------------------------------------------------------------------
      

end
    
  




elseif (steady_state_only == 0)  % Transient runs.
    
    
 % RESET E and fs to values that best fit modern surface and velocity
 set_boundary_conditions_LGM;

 % Grounding line history set in MAIN_Darwin.m
 
    
 for time = 1:N_t_mesh
     
     time2 = time;  % Set this in code, but didn't really need it!
     
   if (time == 1)
       disp('Start in steady state')
   end
   if ( (time == 1) && (min_search_bed == 1) || (min_search_E == 1) || (min_search_fs == 1))
       disp('Start in steady state with minimized parameters')
   end
       disp(' ')
       disp(['Time in calculation = ',int2str(t_P(time))]);
       disp(' ')
       
       
% -------------------------------------------------------------------------    
    calculate_surface_Darwin;
% -------------------------------------------------------------------------

      disp('Calculating for Darwin')

      
% RESET Hatherton from Darwin calculation
disp('  Resetting surface elevation for Hatherton based on Darwin values ... ')
S_in_global_Hatherton = h_w(time2, index_xpos_Hatherton_edges) + B_w(index_xpos_Hatherton_edges);
S_0_in2               = S_in_global_Hatherton;
S_at_GL2(time2)       = S_0_in2;
      
     
% -------------------------------------------------------------------------
    calculate_surface_Hatherton;
      disp('Calculating for Hatherton')
% -------------------------------------------------------------------------
    
    
% Then, setup flux added at 17 km from Hatherton...
flux_add_P(index_xpos_Hatherton_P) = flux_edges_dyn_xt2(time, 1) / W_P2(1);
[ flux_add_w, ...
  flux_add_e ] = get_edge_values_quadratic ...
                                  ( flux_add_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

% -------------------------------------------------------------------------                              
   calculate_surface_Darwin;
      disp('Recalculating for Darwin')
% -------------------------------------------------------------------------

 end   % for loop on time

end  % if statement on steady_state


                 
save save_output.mat




if ((LGM_steady_state ~= 1) && (LGM_transient ~= 1) && ...
    (min_search_E ~= 1) && (min_search_bed ~= 1) && (min_search_fs ~= 1) && (min_search_E_and_fs ~= 1) )


figure(10)   % SURFACE and BED -- Darwin
subplot(2,1,1), plot(x_P/1000, B_P,'r')
hold on
plot(x_P/1000, S_modern, 'linewidth', 2)
plot(x_P/1000, S_P(1,:),'c--')
legend('Bed', 'Measured surface', 'Calculated surface', 'location', 'northwest')
title('Surface and Bed topography: DARWIN')
%xlabel('Distance along flowband (km)')
ylabel('Elevation (m)')
xlim([x_P(1) x_P(end)]/1000)

subplot(2,1,2)
plot(x_P/1000, S_modern-S_P(1,:), 'linewidth', 2)
hold on
plot([x_P(1) x_P(end)]/1000, [0 0],'k--')
title('Modern surface - calculated surface')
xlabel('Distance along flowband (km)')
ylabel('Elevation difference (m)')
xlim([x_P(1) x_P(end)]/1000)



figure(11)   % SURFACE and BED -- Hatherton
subplot(2,1,1), plot(x_P2/1000, B_P2,'r')
hold on
plot(x_P2/1000, S_modern2, 'linewidth', 2)
plot(x_P2/1000, S_P2(1,:),'c--')
legend('Bed', 'Measured surface', 'Calculated surface', 'location', 'northwest')
title('Surface and Bed topography: HATHERTON')
%xlabel('Distance along flowband (km)')
ylabel('Elevation (m)')
xlim([x_P2(1) x_P2(end)]/1000)

subplot(2,1,2)
plot(x_P2/1000, S_modern2-S_P2(1,:), 'linewidth', 2)
hold on
plot([x_P2(1) x_P2(end)]/1000, [0 0],'k--')
title('Modern surface - calculated surface')
xlabel('Distance along flowband (km)')
ylabel('Elevation difference (m)')
xlim([x_P2(1) x_P2(end)]/1000)



figure(30)   % SURFACE VELOCITY 
subplot(2,1,1), plot(x_edges/1000, surf_vel_estimate*(5/4), 'k--', 'linewidth', 2)
hold on
plot(Darwin_measures_centerline_distance/1000, Darwin_measures_flowspeed,'m', 'linewidth', 2)
title('Surface velocity: DARWIN')
ylabel('Velocity (m/yr)')
legend('Calculated velocity',  'MEaSUREs velocity', 'location', 'best')

measures_on_xedges = interp1( Darwin_measures_centerline_distance, Darwin_measures_flowspeed, x_edges, 'linear', 'extrap');

subplot(2,1,2), plot(x_edges/1000, measures_on_xedges - surf_vel_estimate*(5/4), 'b', 'linewidth', 2)
hold on
plot([x_edges(1) x_edges(end)]/1000, [0 0],'k--')
title('Observed - calculated surface velocity')
xlabel('Distance along flowband (km)')
ylabel('Velocity difference (m/yr)')



figure(31)   % SURFACE VELOCITY 
subplot(2,1,1), plot(x_edges2/1000, surf_vel_estimate2*(5/4), 'k--', 'linewidth', 2)
hold on
plot(Hat_measures_centerline_distance/1000, Hat_measures_flowspeed,'m', 'linewidth', 2)
title('Surface velocity: HATHERTON')
ylabel('Velocity (m/yr)')
legend('Calculated velocity',  'MEaSUREs velocity', 'location', 'best')

measures_on_xedges2 = interp1( Hat_measures_centerline_distance, Hat_measures_flowspeed, x_edges2, 'linear', 'extrap');

subplot(2,1,2), plot(x_edges2/1000, measures_on_xedges2 - average_vel_estimate2*(5/4), 'b', 'linewidth', 2)
hold on
plot([x_edges2(1) x_edges2(end)]/1000, [0 0],'k--')
title('Observed - calculated surface velocity')
xlabel('Distance along flowband (km)')
ylabel('Velocity difference (m/yr)')



figure(40)
plot(x_edges, flux_edges_dyn_xt(1,:))
hold on
plot(x_edges, flux_edges_kin_xt(1,:),'c')
plot(x_edges, ones(size(x_edges)),'k--')
title('DARWIN')


figure(41)
plot(x_edges2, flux_edges_dyn_xt2(1,:))
hold on
plot(x_edges2, flux_edges_kin_xt2(1,:),'c')
plot(x_edges2, ones(size(x_edges2)),'k--')
title('HATHERTON')



end


% % Can this really be used to evaluate relative importance of dynamics vs.
% % mass balance ... holds for ice mass frozen to the bed.
% BDM_velocity      = interp1(measures_centerline_distance, measures_flowspeed, [x_P]);
% F_check = ( (S_modern - B_P)./b_dot_P(1,:)) .* (BDM_velocity ./ x_P);







