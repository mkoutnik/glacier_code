function   [ h_P_t, h_w_t, h_e_t, ...
             S_P_t, S_w_t, S_e_t, ...   
             flux_edges_dyn_t   ] = explicit_solver( h_P, S_P, ...
                                       x_P, x_w, x_e, ...
                                       dx_P, dx_w, dx_e, ...
                                       B_P, B_w, B_e, W_P, W_w, W_e,  ...
                                       slip_P, slip_w, slip_e, ...
                                       A_eff_edges, b_dot_P, ...
                                       b_dot_edges_SS, b_dot_edges, ...
                                       Q_out_L_SS, Q_out_R_SS, ...
                                       Q_external_L, Q_external_R, ...
                                       S_at_GL, ...
                                       t_P, dt_P, i_time)


%--------------------------------------------------------------
%
% Michelle Koutnik (mkoutnik@ess.washington.edu)
% adapting and expanding pieces of code from Ed Waddington
% last updated, June 2008
%
% Explicit solution to continuity equation for ice-sheet surface evolution

% Detailed comments at the END of this file.

%---------------------------------------------------------------



% '_t' values are those that are unknown, and being solved for at the 
% time = time; known values are at time-1 and have no underscore

% all values sent in to the implicit solver are the known
% values at time-1, solver will find the unknown values at time=time 




% known values
% ============
  [ h_w, h_e ] = get_edge_values_quadratic( h_P, x_P, x_w, x_e, ...
                                            dx_P, dx_w, dx_e );
                     
    
  [ dS_dx_w, dS_dx_e ] = get_gradient_values( S_P, x_P, dx_P );
  
  
  
% Boundary condition
% ==================                     
                      
% dynamic flux. 
% ------------- 

Q_w = calc_flux_dyn( x_w, h_w, dS_dx_e, slip_w, W_w, ...
                            A_eff_edges(1:end-1) );  

Q_e = calc_flux_dyn( x_e, h_e, dS_dx_e, slip_e, W_e, ...
                            A_eff_edges(2:end) );  
                         

% Q_in_use = Q_out_L_SS - Q_external_L(1);  
% Q_w = [ Q_in_use Q_w(2:end) ];                       
       
 
 b_dot_w = b_dot_edges(1, 1:end-1);
 b_dot_e = b_dot_edges(1, 2:end);
 
       
 W_bar = (W_e + W_w)/2;
 
 
 
 % dS/dt, thickness change for calculating S_P at future timestep
 % use the flux that is known for the current timestep
 % and the accumulation rate at future timestep
 % ---------------------------------------------------
     S_dot_i_time = - ( (Q_e - Q_w ) ./ dx_P  ...
                    - ( 2* (b_dot_w .* W_w + b_dot_e .* W_e  )   ...
                    + (b_dot_w .* W_e + b_dot_e .* W_w  ) ) / 6) ./W_bar;
              
               
               
 % S, ice-surface elevation at future timestep
 % -------------------------------------------
        S_P_t = (S_P) + (dt_P * S_dot_i_time);
 
        h_P_t = (S_P_t - B_P);

        if (sum( isnan(S_P_t) > 1))
            check_here = 1;
        end


   [ h_w_t, h_e_t ] = get_edge_values_quadratic( h_P_t, x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );

   
                                             
% recalculate flux
% ================
                                             
   [ dS_dx_w_t, dS_dx_e_t ] = get_gradient_values( S_P_t, x_P, dx_P );
                         
   [ S_w_t, S_e_t ] = get_edge_values_quadratic( S_P_t, x_P, x_w, x_e, ...
                                                 dx_P, dx_w, dx_e );                                      
   
   Q_w_t = calc_flux_dyn( x_w, h_w_t, dS_dx_e_t, slip_w, W_w, ...
                            A_eff_edges(1:end-1) );  

   Q_e_t = calc_flux_dyn( x_e, h_e_t, dS_dx_e_t, slip_e, W_e, ...
                            A_eff_edges(2:end) );  
                                    
                        
   flux_edges_dyn_t = [ Q_w_t(1) Q_e_t ]; 
    
   
   
   
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETAILED COMMENTS:
%  explicit scheme - really easy!  (if it works ...)
%  finds transient ice-sheet surface elevation S(x,t)
%
%  S_0(T_0)   - initial condition
%  b_dot      - accumulation pattern interpolated onto (x_P, t_P) mesh
%  Q_in       - flux at left edge of control volume at left
%                margin of calculation grid
%  B_P(x)     - bed elevation at grid points x_P
%  W_P(x)     - width of control volume at x_P
%  W_w, W_e   - widths at upstream and downstream edges of control volumes
%
%  (x_P, t_P) - calculation mesh for surface S(x,t) 
%  x_w, x_e   - upstream and downstream edges of control volumes
%  dx_P       - lengths of control volumes (distance between centers of
%               intervals to upstream and downstream x_P points
%  dx_w, dx_e - distances from x_P to upstream and downstream x_P points
%
%-----------------------------------------------------------------------


