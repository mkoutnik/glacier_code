function [ x_P, x_w, x_e,    ...
           dx_P, dx_w, dx_e, ...
           t_P, dt_P,        ...
           z_hat, dz_hat,    ...
           x_edges  ] =      ...
                             load_mesh( x_nodes, t_nodes, solver_method )


% solver_method is whether the code is going to be setup
% for an explicit or implicit solution, a different grid
% resolution is required in these two cases                         
                         

%--------------------------------------------------------------
%
%  Define space and time mesh for numerical solution S(x,t)
%  Create image points outside x_nodes at both ends where 
%  fluxes will be calculated or assigned
%  Surface S(x,t) is found at times t_P and positions x_P

%---------------------------------------------------------------


global min_search_E min_search_fs min_search_bed min_search_E_and_fs
global lower_resolution



if ( (lower_resolution == 1) || (min_search_E == 1) || (min_search_bed == 1) || (min_search_fs == 1) || (min_search_E_and_fs == 1))
    
    N_t = 201;
    N_x = 76;
    
else
    
    N_t = 201;
    N_x = 76; % 151; % 76; 
    
end
    
    

  t_span = t_nodes(end) - t_nodes(1);
  dt = t_span/( N_t - 1 );
  
  t_P = ( t_nodes(1): dt: t_nodes(end) )';
  
  dt_P_temp = diff(t_P)';
  dt_P = [ dt_P_temp dt_P_temp(end) ]';
  
  
  
  x_nodes_span = x_nodes(end) - x_nodes(1);
  dx = x_nodes_span/(N_x - 1);
 
  
  
%  It is better to define the x_P points first, then to infer the
%  control-volume edges. Defining edges first can lead to oscillations 
%  in the centering of x_P in the control volumes
%  Set the first 2 points at each end
%  ===================================================================

% NOTE! There are many calls to interp1q, so use uniformly spaced x_P

      x_P = zeros( 1, N_x-1);
      x_P(1:2)   = x_nodes(1)       + (dx/2) * [1 3];
      x_P(end-1:end) = x_nodes(end) - (dx/2) * [3 1];
   
  
%  Here is the place to define irregular x-spacing (if needed)


% For now, use uniform spacing:
% =============================
    x_P(3:end-2) = [x_P(2)+dx: dx: x_P(end-1) - dx/2 ]; 
%

%   Distances to upstream and downstream P points are given by
%   dx_w and dx_e
      dx_w = diff( x_P );
  
  % add first dx_w = dx for x_P(1)
      dx_w = [ dx_w(1) dx_w ];  
      
      dx_e = diff( x_P );
      
  % add last dx_e = dx for x_P(end)
      dx_e = [ dx_e dx_e(end) ];

   
%   Fluxes Q(x,t) and width W(x) are defined at 
%   upstream and downstream edges x_w and x_e of control volumes 
      x_w = x_P - dx_w/2; 
      x_e = x_P + dx_e/2; 

%  Lengths of control volumes are dx_S
      dx_P = x_e - x_w;


 % group values at edges into one vector
 % =====================================
     x_edges = [ x_w(1) x_e ];     
      
      

% setup vertical grid, z_hat = (z-B(x))/(S(x,t)-B(x))
% nondimensional z grid, z=S=1, z=B=0
% ===================================================
    dz_hat = 0.005; % 0.002 gives the same result
    z_hat  = [ 1: -dz_hat: 0 ];     

    
    
