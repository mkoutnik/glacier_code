function  QQ = calc_flux_dyn( x, h, dS_dx, E, fs, W, A_T, deformation_only, deformation_plus_sliding, sliding_only )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  finds fluxes based on thickness h(x), slope dS_dx, 
%    basal slip r atio, and temperature
%  W(x) is flow-band width
%  A_T is softness parameter in Glen flow law
%---------------------------------------------------------


%  get depth-averaged velocities
     u_bar = get_u_bar( x, h, dS_dx, E, fs, A_T, deformation_only, deformation_plus_sliding, sliding_only );
%

%   form flux
     QQ = h .* u_bar .* W;
%

%
%
%

