function [ T_surf ] = load_T_surf (t_P)



global N_x_mesh  N_t_mesh


% estimate surface temperature in Kelvin!


  T_surf_0 = zeros(N_t_mesh, N_x_mesh);

  
  surface_T = 273.15 - 25;  % from "Elaine" weather station; Reusch and Alley (2004)                                               
  T_surf    = surface_T + T_surf_0;


  
  