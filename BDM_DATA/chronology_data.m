function S_at_GL = chronology_data (t_P) 

global S_in_global


S_0_in = S_in_global;


% Use Perry's history:
 % --------------------
   load mt_hope_exposure_ages.mat  % exposure_age_yrs_Hope elevation_m_asl_Hope ...
                                   % internal_uncertainty_yrs_Hope external_uncertainty_yrs_Hope
   relative_elevation_Hope = S_0_in + (elevation_m_asl_Hope - min(elevation_m_asl_Hope)); 
   S_at_GL                 = interp1(-exposure_age_yrs_Hope, relative_elevation_Hope, t_P, 'linear', 'extrap');
   S_at_GL(end)            = S_0_in; % reset to modern elevation
 