function [ val_w, val_e ] = get_edge_values_quadratic ...
                               ( val_P, x_P, x_w, x_e, dx_P, dx_w, dx_e )



% -----------------------------------------------------------------------
% from Ed Waddington.

%   From values at val_P at interior points x_P
%   interpolate to get ice thickness at edges of control volumes
%   Extrapolate values on ends for points beyond mesh.

%   Since I use a linear extrapolation for slopes at the 
%   control-volume faces that are off the x_P grid, I should
%   use a quadratic extrapolation for the ice thickness.

% -----------------------------------------------------------------------


%  get values at upstream edges of control volumes
%  --------------------------------------------

%  first extrapolate to imaginary downstream edge of first c.v.


     c_1 = dx_w(1)/(2*dx_e(1) );
     c_2 = dx_w(1)/(4*dx_P(2) );
     c_3 = dx_w(1)/(2*dx_e(2) );
     c_4 = dx_e(1)/(2*dx_P(2) );

 
     val_w_0 = val_P(1) * ( 1 + c_1 * (1 + c_2 + c_4 ) ) ...
              - val_P(2) * ( c_1 * (1 + c_2 +c_4 ) + c_3 * ( c_2 + c_4 ) ) ...
              + val_P(3) * ( c_3 * ( c_2 + c_4 ) );

   % val_w = [  val_w_0  qinterp1( x_P, val_P, x_w(2:end), 0 )' ];
     val_w = [  val_w_0  interp1( x_P, val_P, x_w(2:end) ) ];

     
     
%  get points at downstream edges control volumes
%  --------------------------------------------                

%  first extrapolate to imaginary downstream edge of last c.v.


     c_1 = dx_e(end)/(2*dx_w(end) );
     c_2 = dx_e(end)/(4*dx_P(end-1) );
     c_3 = dx_e(end)/(2*dx_w(end-1) );
     c_4 = dx_w(end)/(2*dx_P(end-1) );
                
  
     val_e_end = val_P(end) * ( 1 + c_1 * (1 + c_2 + c_4) ) ...
                 - val_P(end-1) * ( c_1 * (1 + c_2 + c_4 ) + ( c_2 + c_4 ) * c_3 ) ...
                 + val_P(end-2) * ( ( c_2 + c_4 )* c_3 );

    % val_e = [ qinterp1( x_P, val_P, x_e(1:end-1), 0 )'  val_e_end ];
      val_e = [ interp1( x_P, val_P, x_e(1:end-1) )  val_e_end ];
 
     
     
     
     