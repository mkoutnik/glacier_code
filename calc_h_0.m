 function [ S_0, h_0, ...
            dS_dx_P, ...
            dS_dx_edges ] = calc_h_0( x_P, x_w, x_e, dx_P, ...
                                      B_P, B_w, B_e, W_P, W_w, W_e, ...
                                      b_dot_edges, E_w, E_e, fs_w, fs_e, ...
                                      Q_0_in, S_0_in, A_T_w, A_T_e, ...
                                      flux_add_w, flux_add_e, ...
                                      deformation_only, deformation_plus_sliding, sliding_only)

     
%--------------------------------------------------------------------------
% from Ed Waddington
%
% This code calculates initial condition for transient surface calculation
% Detailed comments at the END of this file.
%
%--------------------------------------------------------------------------


% edge values of accumulation rate
% ================================
  b_dot_w = b_dot_edges(1:end-1);
  b_dot_e = b_dot_edges(2:end);

       
%  get length of x mesh      
       N_x = length(x_P);

%  set maximum number of allowable iterations with Newton's method
       itmax = 5000;
        
%  initialize residuals and thickness-change matrices for Newton's method 
       r_save  = NaN( itmax, N_x );
       dS_save = NaN( itmax, N_x );

%  set tolerance test for final residuals
       r_limit = 1e-8;
       
%  set fractional thickness perturbation when evaluating 
%  derivatives dr/d S_e
       h_frac = 1e-5;
   
   
%  Initialize surface-profile vectors
      S_w     = zeros( 1, N_x );
      S_e     = S_w;
      S_0     = S_w;
      slope_w = zeros( 1, N_x+1 );


%   set elevation to S_0_in at x_w(1)
      S_w(1)  = S_0_in;

%  use known thickness at x_w(1) to find slope at x_w(1)
      slope_w(1) =  dSurf_0_dx( x_w(1), x_P, x_w, x_e, dx_P, ...
                               S_w(1), B_w, B_e, W_w, W_e, ...
                               E_w, E_e, fs_w, fs_e, b_dot_w, b_dot_e, ...
                               Q_0_in, A_T_w, A_T_e, flux_add_w, flux_add_e, ...
                               deformation_only, deformation_plus_sliding, sliding_only);                        
%

%   Now get surface profile S_0( x_P ) consistent with numerical mesh 
%   -----------------------------------------------------------------

  for i = 1:N_x
      
    %  extrapolate slope at x_w(i) to find S_0(i) at x_P(i) 
        dx  = (x_P(i) - x_w(i)); 
        S_0(i) = S_w(i) + slope_w(i) * dx;
    %
 
    %  Find slope at S_0(i) and extrapolate to get first estimate for
    %  surface S_w(i) at next control-volume edge
    %  (a) get slope
        slope_i = dSurf_0_dx( x_P(i), x_P, x_w, x_e, dx_P, ...
                               S_0(i), B_w, B_e, W_w, W_e, ...
                               E_w, E_e, fs_w, fs_e, b_dot_w, b_dot_e, ...
                               Q_0_in, A_T_e, A_T_w, flux_add_w, flux_add_e, ...
                               deformation_only, deformation_plus_sliding, sliding_only);                        

    
    %  (b) extrapolate surface to x_e(i)
        dx = (x_e(i) - x_P(i));
        S_e_trial = S_0(i) + slope_i * dx;
  
    
  %  Use Newton's method to find "correct" thickness at x_e(i)
  %  The correct thickness produces a steady-state slope at x_e(i)
  %  that projects back to S_0(i) with small residual "r" 
  %     r = S_0(i) - ( S_e_trial - slope_e * dx )
  %  
  %  Numerically estimate dr/d S_e.  Correction dS to surface at x_e(i)
  %  that would eliminate residual (if system was linear) is given by
  %     dr/d S_e * dS = -r
  %  and S_e <-- S_e + dr
  
  %  reset counter for iterations
        itno = 1;
        r_a = 99999999;
  
       
  %  loop to find corrections to S_e(i)
       while ( (itno <= itmax) & (abs(r_a) > r_limit) )
  
           
     %  find trial slope at x_e(i)
          slope_e_a =  dSurf_0_dx( x_e(i), x_P, x_w, x_e, dx_P, ...
                                   S_e_trial, B_w, B_e, W_w, W_e, ...
                                   E_w, E_e, fs_w, fs_e, b_dot_w, b_dot_e, ...
                                   Q_0_in, A_T_e, A_T_w, flux_add_w, flux_add_e, ...
                                   deformation_only, deformation_plus_sliding, sliding_only);                        

  
     %  project back to x_P(i)
          r_a = S_0(i) - ( S_e_trial - dx * slope_e_a );
  
     %  perturb thickness at x_e(i)
          dh = h_frac * ( S_e_trial - B_e(i) );
  
  
     %  find trial slope at x_e(i) with perturbed thickness
          slope_e_b =  dSurf_0_dx( x_e(i), x_P, x_w, x_e, dx_P, ...
                               (S_e_trial+dh), B_w, B_e, W_w, W_e, ...
                               E_w, E_e, fs_w, fs_e, b_dot_w, b_dot_e, ...
                               Q_0_in, A_T_e, A_T_w, flux_add_w, flux_add_e, ...
                               deformation_only, deformation_plus_sliding, sliding_only);                        


     %  project back to x_P(i)
          r_b = S_0(i) - ( (S_e_trial + dh) - dx * slope_e_b );
     
  
     %  estimate derivative dr/dS_e
          dr_dS = (r_b - r_a)/dh;
     

     %  find correction dS to S_e(i)
          dS = -r_a/dr_dS;
     
     %  get new estimate of S_e(i)
          S_e_trial = S_e_trial + dS;
     %
     %  save current residual and current dS
          r_save(itno, i )  = r_a;
          dS_save(itno, i ) = dS;
       
       
     % get ready for next iteration
          S_e(i)   = S_e_trial;
          S_w(i+1) = S_e_trial;
          slope_w(i+1) = slope_e_a;
          
          itno = itno+1;
       
        if (itno == itmax)
            disp(' Hit maximum iterations in calc_h_0.m')
        end
          
      end   % while( (itno < itmax) & ...
       
      
      
  end   %  for i = 1:N_x ...



% ice thickness
% =============
  h_0 = S_0 - B_P;
  
  
 
% surface slope
% =============
 [ dS_dx_w, dS_dx_e ] = get_gradient_values( S_0, x_P, dx_P );
 
   dS_dx_edges = [dS_dx_w(1) dS_dx_e];
   dS_dx_P     = qinterp1([x_w(1) x_e], dS_dx_edges, x_P)';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETAILED COMMENTS:

% from Ed Waddington

%  Set up initial (steady-state) condition S_0(x) for ice-surface height.
%  Routine uses a 2-step algorithm that is consistent with numerical
%  scheme, in the sense that it should produce a surface profile that will
%  be steady in the numerical scheme with the specified mesh.
%
%  Surface elevation must be specified at one point.
%  Boundary condition S_0_in is applied at "upstream" (left) boundary,
%  x_w(1), which is a flux point rather than a thickness point.
%
%  When thickness is known at a flux point x_w(i), we can calculate the
%  slope there and project to find surface S_0(i) at next "downstream"
%  thickness point x_P(i).  (This linearity between thickness points is
%  assumed assumption in numerical model.)
%
%  Then we need to find thickness S_e(i) at next "downstream" control-volume
%  edge x_e(i) such that the slope there, when projected back, exactly
%  reproduces S_0(i) at x_P(i), which again is the assumption of linearity
%  between thickness points in the numerical scheme.
%  Routine uses Newton's method to find best S_e(i). 
%
%
%
%   In this first simple steady-state model, we use Shallow Ice
%   Approximation to integrate dS/dx with flux derived kinematically 
%   from accumulation b_dot(x,t) and width W(x)
%
%   This routine could be generalized so that thickness could be
%   specified at any x, rather than always at "upstream" boundary x_w(1).
%   But not to bother now.

