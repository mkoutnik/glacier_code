   function  X = find_X_intervals( x_flux, x_mesh )
   
   
   
% -------------------------------------------------------------------------
% from Ed Waddington

%   set up nondimensional array X(i,j)
%   X(i,j) = (x_ij - x_i) / (x_i+1 - x_1) = (x_ij- x_i) / dx_i 
%   column j gives the  nondimensional value X of the point at which the
%   integration ends within each x_mesh-defined interval [ x_i  x_i+1 ], when
%         Q(x_flux(j) ) = integral_{x_1}^{x_flux(j)}  [ b_dot(x) W(x) ]  dx
%
%   x_i are the x_mesh values x_mesh(i) at which b_dot and W are specified
%   x_flux(j) are the points at which we want to find Q( x_flux(j) )
%   x_ij are the points at which integration ends in
%   interval i =  [ x_i  x_i+1 ] when finding  Q(x_j).    
%   x_ij = x_i+1  when x_i+1 < x_flux_j
%             = x_j     when x_i < x_flux_j  < x_i+1
%             = x_i     when x_flux_j < x_i
%  But it is easier to deal with nondimensionalized positions within
%  each x_nodes interval   [ x_i  x_i+1 ] : 
%  ------------------------------------------------------------------------


%  get number of mesh points (where b_dot and W are defined)
      N_mesh = length( x_mesh );
      

%  get lengths of intervals dx_i between x_mesh(i) and x_mesh(i+1)
%  -----------------------------------------------------------------
      dx_temp = diff( x_mesh );

   %  add an extra interval to accommodate requests for flux at x_mesh(end)
   %  ----------------------------------------------------------------------
          dx = [ dx_temp dx_temp(end) ];


%  find which x_mesh interval grid_boxes(j) contains each x_flux(j)
%  ----------------------------------------------------------------
  
 %   grid_boxes = where_is_x( x_flux, x_mesh );


% put this function inline for speed:
  
%  clunky but simple search
%  -------------------------
       N_x_flux = length( x_flux );
       grid_boxes = NaN(1,N_x_flux);
       
      for j=1:N_x_flux
         index = find( x_mesh <= x_flux(j) );
         grid_boxes(j) = index(end);
      end  %  for j=1:N_x ...

    
    
    
    
%  find nondimensional position X within its grid box for each x_flux(j)
%  ---------------------------------------------------------------------
      X_j =  ( x_flux - x_mesh( grid_boxes) ) ./ dx( grid_boxes );


% ------------------
%  STUFF X(i,j)
% ------------------

   %  initialize X with zeros   (no contribution to flux)
   %  ---------------------------------------------------
         X = zeros( N_mesh, length(x_flux) );


   %  vector "index_0" indexes the locations of the first interval at x_in (top row) 
   %  vector "index_1" indexes the locations of the partial intervals X_j
   %  -------------------------------------------------------------------
       index_0 = ( cumsum( N_mesh * ones( size(x_flux) ) ) - N_mesh)+1;
       index_1 = (index_0 -1) + grid_boxes;


   %  now stuff X(i,j) = 1 into rows grid_boxes(j) and above in column j
   %  i.e. the intervals that contribute fully to Q(x_j)
   %  ---------------------------------------------------- 
        for j = 1:length(x_flux)
             index = index_0(j):index_1(j);
             X(index) = ones( size(index) );
         end     % for j = 1:length(x_flux) ...


   %  put the partial integration lengths X_j into the locations index_1(j),
   %   which correspond to  X( grid_boxes(j), j )
   %  --------------------------------------------------------
          X(index_1) = X_j;

%
%
%
%

