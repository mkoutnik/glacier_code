% run_min_search_GL_at_LGM


% Load Hatherton LGM data:
load LGM_values_use.mat

  disp (' ')
  disp('Assuming that LGM maximum at the same time at all locations of data... ')
  disp(' ')  
% DAN_ages_use DAN_elevations_use DAN_position 
% MVfloor_ages_use MVfloor_elevations_use MVfloor_position
% LWC14_ages_use LWC14_elevations_use LWC14_position
    

LGM_elev_vec = round(S_modern(1)) + [0:20:1000];
N_values     = length( LGM_elev_vec );

S_0_in_orig  = S_0_in;

RMS_mismatch = NaN * ones( 1, N_values );



 for i_values = 1:N_values   % loop over all x positions -- downstream
    
disp(' ')     
disp(['Evaluating LGM estimate of ',int2str(LGM_elev_vec(i_values)), ' meters, which is ', int2str(i_values), ' of total estimates = ',int2str(N_values)]);

     
% Value for Darwin GL
    S_0_in = LGM_elev_vec (i_values);
    S_at_GL(time) = S_0_in;
                         
 
% -------------------------------------------------------------------------    
    calculate_surface_Darwin;
% -------------------------------------------------------------------------

      disp('Calculating for Darwin')

      
% RESET Hatherton from Darwin calculation
disp('  Resetting surface elevation for Hatherton based on Darwin values ... ')
%S_in_global_Hatherton = h_w(1, index_xpos_Hatherton_edges) + B_w(1, index_xpos_Hatherton_edges) - value_to_zero_bed + value_to_zero_bed2;
S_in_global_Hatherton = h_w(1, index_xpos_Hatherton_edges) + B_w(1, index_xpos_Hatherton_edges);
S_0_in2               = S_in_global_Hatherton;
S_at_GL2              = ones(size(t_P2)) * S_0_in2;  
      
     
% -------------------------------------------------------------------------
    calculate_surface_Hatherton;
      disp('Calculating for Hatherton')
% -----------------------------------------
 
 
 
 % Residual
 % ========
    std_dev = 1; 

 % Compare Hatherton values to data:
    DAN_mean_age     = round(mean(DAN_ages_use));
    MVfloor_mean_age = round(mean(MVfloor_ages_use));
    LWC14_mean_age   = round(mean(LWC14_ages_use)); 
 
   DAN_mean_elevation = mean(DAN_elevations_use);
   DAN_calc_elevation = interp1( x_P2, S_P2(time,:), DAN_position ); 
   residual_DAN = sqrt( mean( ( ( abs(DAN_mean_elevation) - abs(DAN_calc_elevation) )/std_dev).^2) );
 
   MVfloor_mean_elevation = mean(MVfloor_elevations_use);
   MVfloor_calc_elevation = interp1( x_P2, S_P2(time,:), MVfloor_position ); 
   residual_MVfloor = sqrt( mean( ( ( abs(MVfloor_mean_elevation) - abs(MVfloor_calc_elevation) )/std_dev).^2) );
   
   LWC14_mean_elevation = mean(LWC14_elevations_use);
   LWC14_calc_elevation = interp1( x_P2, S_P2(time,:), LWC14_position ); 
   residual_LWC14 = sqrt( mean( ( ( abs(LWC14_mean_elevation) - abs(LWC14_calc_elevation) )/std_dev).^2) );
   
   
   
   residual = residual_DAN + residual_MVfloor + residual_LWC14;
    
  
   RMS_mismatch(i_values) = residual;
   
                                   
 end  % loop over E values
 
 
   index     = find( RMS_mismatch == min(min(RMS_mismatch) ) );
   
%    if (length(index)>1)  % Multiple values of E give same mismatch
%        disp('Multiple values of E give same mismatch value -- set E=1')
%        index = find(E_vec == 1);
%    end
   
   RMS_best  = RMS_mismatch(index);
   best_LGM_elev  = LGM_elev_vec(index)

   

 disp(' Running again with best LGM elevation at Darwin GL ... ')  
   
     
% Value for Darwin GL
    S_0_in = best_LGM_elev;
    S_at_GL(time) = S_0_in;
                         
 
% -------------------------------------------------------------------------    
    calculate_surface_Darwin;
% -------------------------------------------------------------------------

      disp('Calculating for Darwin')

      
% RESET Hatherton from Darwin calculation
disp('  Resetting surface elevation for Hatherton based on Darwin values ... ')
%S_in_global_Hatherton = h_w(1, index_xpos_Hatherton_edges) + B_w(1, index_xpos_Hatherton_edges) - value_to_zero_bed + value_to_zero_bed2;
S_in_global_Hatherton = h_w(1, index_xpos_Hatherton_edges) + B_w(1, index_xpos_Hatherton_edges);
S_0_in2               = S_in_global_Hatherton;
S_at_GL2              = ones(size(t_P2)) * S_0_in2;  
      
     
% -------------------------------------------------------------------------
    calculate_surface_Hatherton;
      disp('Calculating for Hatherton')
% -----------------------------------------
 

% Then, setup flux added at set km from Hatherton...
flux_add_P(index_xpos_Hatherton_P) = flux_edges_dyn_xt2(time, 1) / W_P2(1);
[ flux_add_w, ...
  flux_add_e ] = get_edge_values_quadratic ...
                                  ( flux_add_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

% -------------------------------------------------------------------------                              
   calculate_surface_Darwin;
      disp('Recalculating for Darwin')
% -------------------------------------------------------------------------
      




% FIGURES

% mismatch vs. LGM elevation values
figure(100)
plot(LGM_elev_vec, RMS_mismatch, 'k', 'linewidth', 2)
set (gca, 'fontsize', 14)
xlabel('LGM elevation at Darwin Glacier (m)', 'fontsize', 16, 'fontweight', 'bold')
ylabel({'RMS mismatch between model'; 'and data along Hatherton Glacier'}, 'fontsize', 16, 'fontweight', 'bold')



% Elevation profiles along Hatherton vs. data locations and elevation
figure(101)   % SURFACE and BED -- Hatherton
subplot(2,1,1), plot(x_P2/1000, B_P2,'r')
hold on
plot(x_P2/1000, S_modern2, 'linewidth', 2)
plot(x_P2/1000, S_P2(1,:),'c--')
plot(DAN_position/1000, DAN_mean_elevation, 'g*')
plot(MVfloor_position/1000, MVfloor_mean_elevation, 'm*')
plot(LWC14_position/1000, LWC14_mean_elevation, 'k*')
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






