function u_bar = get_u_bar( x, h, dS_dx, E, fs, A_T, deformation_only, deformation_plus_sliding, sliding_only )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  form depth-averaged velocity at positions x, using
%   thickness h(x), slope dS_dx(x), basal slip ratio,
%   and ice softness A_T(x) at effective temperature
%----------------------------------------------------


     global rho_ice  g  n

%      global deformation_only deformation_plus_sliding sliding_only ...
%        deformation_sliding_lateraldrag deformation_sliding_longstress ...
%        deformation_sliding_lateraldrag_longstress
%    
   
     
%  use Shallow Ice Approximation with effective temperature

if (deformation_only == 1)
     
    u_bar = ( -sign(dS_dx) .* ( (2*E.*A_T)./(n+2) ) .* ...
               (rho_ice * g * abs(dS_dx) ).^n ) .* h.^(n+1);

           
           
elseif (deformation_plus_sliding == 1)

    deformation_term = ( -sign(dS_dx) .* ( (2*E.*A_T)./(n+2) ) .* ...
               (rho_ice * g * abs(dS_dx) ).^n ) .* h.^(n+1);
    
    sliding_term = ( -sign(dS_dx) .* (fs .* (rho_ice * g * abs(dS_dx) ).^n ) .* h.^(n-1) );
    
    u_bar = (deformation_term + sliding_term);
    
    
    
elseif (sliding_only == 1)
    
    deformation_term = 0;
    
    sliding_term = ( -sign(dS_dx) .* (fs .* (rho_ice * g * abs(dS_dx) ).^n ) .* h.^(n-1) );
    
    u_bar = (deformation_term + sliding_term);
    
end
%
%
%

