function dSurf_0_dx = dSurf_0_dx_flux_add( x, x_P, x_w, x_e, dx_P, S, ...
                                  B_w, B_e, W_w, W_e, slip_w, ...
                                  slip_e, b_dot_w, b_dot_e, ...
                                  Q_in_0, A_T_w, A_T_e, flux_add )
                        

% -------------------------------------------------------------------------
% from Ed Waddington
%
%   Evaluates derivative dS/dx at x, given S(x), when integrating dS/dx
%   to find steady-state surface profile S_0(x)
%   with accumulation pattern b_dot(x) and incoming flux Q_in_0
%     x can be a row vector.
%
%   Flux is calculated kinematically
%     Q(x,t) = Q_in_0 + int_{x_in}^x [b_dot(x') - S_dot(x')] W(x') dx'
%
%   with S_dot(x) = 0  for steady state
%   Then slope dS/dx is calculated from SIA
%    dS/dx(x) = (1/ rho g ) ...
%          [ [ (n+2) Q(x)/( 2 A(T) )] 1/(S(x)-B(x) )^(n+2) ]^(1/n)
%
% -------------------------------------------------------------------------
  
global rho_ice  g  n
   

     
%  find bed elevation, flow-band width, and slip ratio at x
%    (Recall that bed is defined primarily at c-v interior points
%    while width and slip are defined at control-volume edges)
%  ------------------------------------------------------------
     Bed    = interp1( [ x_w(1)  x_e ], [ B_w(1)  B_e ], x ); 
     Width  = interp1( [ x_w(1)  x_e ], [ W_w(1)  W_e ], x );
     slip_x = interp1( [ x_w(1)  x_e ], [slip_w(1) slip_e], x );
     
     E_x = slip_x;
     
     A_T_x = interp1( [ x_w(1)  x_e ], [ A_T_w(1) A_T_e ], x ); 

%  set dS/dt = 0  because we are finding a SS profile
      S_dot_0 = zeros( size([x_w(1) x_e]) );
     
      
%  find the flux at points x
%  -------------------------
     Q_x = calc_flux_kin( x, [x_w(1) x_e], [dx_P dx_P(end)], [W_w(1) W_e], ...
                          S_dot_0, [b_dot_w(1) b_dot_e], Q_in_0 );
%
     Q_x = Q_x + flux_add;

%  find the surface slope
%  ----------------------
     dSurf_0_dx = - sign(Q_x) .* (1/(rho_ice*g) ) .* ...
             ( ( (n+2)./(2*E_x.*A_T_x) ) .* ( abs(Q_x) ./ Width ) ...
                               ./ ((S - Bed ).^(n+2) ) ).^(1/n);                  

%
%
%

