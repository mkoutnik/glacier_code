function [ dvalue_dx_w, dvalue_dx_e ] = ...
                                    get_gradient_values( value, x_P, dx_P )


% -------------------------------------------------------------------------
% from Ed Waddington 

%   find slopes at upstream, downstream, and center points
% -------------------------------------------------------------------------


x_w     = x_P - (dx_P/2);
x_e     = x_P + (dx_P/2);
x_edges = [x_w(1) x_e];


%  form slope vector dvalue_dx_P
  %   x_P_mat = repmat(x_P, size(value,1),1);
      x_P_mat = x_P( ones(size(value,1),1), :); 
  
     dvalue_dx = diff(value) ./ diff(x_P_mat);
       
 
     
%  make 2 vectors, one for upstream edges, one for downstream edges
       
%  upstream edges of control volumes
   %  first extrapolate linear slope variation to upstream edge
   %  of first c.v.
   % ==========================================================
      dvalue_dx_0 = dvalue_dx(1) * (1 + dx_P(1)/dx_P(2) ) ...
                    - dvalue_dx(2) * (dx_P(1)/dx_P(2) );
     
      dvalue_dx_w = [  dvalue_dx_0   dvalue_dx  ];

      
%  downstream edges of control volumes
   %  first extrapolate to downstream edge of last c.v.
   %  (remember that dS_dx has one element fewer than S and dx_P)
   % ============================================================
      dvalue_dx_end = dvalue_dx(end) * (1 + dx_P(end)/dx_P(end-1) ) ...
                    - dvalue_dx(end-1) * ( dx_P(end)/dx_P(end-1) );

      dvalue_dx_e = [ dvalue_dx  dvalue_dx_end ];


    




