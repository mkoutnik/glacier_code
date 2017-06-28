% script size_nodes_mesh.m


%--------------------------------------------------------------------------
%
% Michelle Koutnik (mkoutnik@ess.washington.edu)
% last updated, January 2009
%
% This script establishes and loads in global variables that are used
% throughout the code.

% sets up the size of the nodal and mesh-grid values
% these are GLOBAL values


%--------------------------------------------------------------------------

global N_t_nodes  N_x_nodes  N_t_mesh  N_x_mesh   N_z


% dt_mesh    = timestep
% N_t_nodes  = number of time points
% N_x_nodes  = number of spatial points
% N_t_mesh   = number of timesteps
% N_x_mesh   = number of spatial points
% N_data     = number of data
% N_z        = number of vertical points



% number of x_nodes and t_nodes
     N_t_nodes = length (t_nodes);
     N_x_nodes = length (x_nodes);

     
%  N_t is number of time steps
     N_t_mesh = length(t_P);

     
%  N_x is number of control volumes
     N_x_mesh = length(x_P);

     
% N_z is the number of z values in the grid
     N_z = length(z_hat);
     
     
     
     