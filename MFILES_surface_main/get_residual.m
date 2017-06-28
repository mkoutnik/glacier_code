function [ residual ] = get_residual( h_P, h_P_t, h_P_t_guess, ...
                                      x_P, x_w, x_e, dx_P, dt_P, ...
                                      Q_w, Q_e, Q_w_t, Q_e_t, ...
                                      K_w_t, K_e_t, ...
                                      W_P, W_w, W_e, b_dot_P, ...
                                      Q_out_L, Q_out_L_t, ...
                                      Q_out_R, Q_out_R_t, ...
                                      full_domain_flag  )
                                      

% ----------------------------------------------------------------------- %
% Michelle Koutnik
% last updated: January 2009
%
% Detailed comments at the END of this file.
%
% ----------------------------------------------------------------------- %
                                 
           

global theta 

 


% residual of continuity equation
% with new estimate of h_P_t
% ==========================


  residual = ((h_P_t - h_P) .* (dx_P./dt_P)) + ...
             ( (1./W_P) .* ( (theta* (Q_e_t - Q_w_t)) + ...
                            ((1-theta)* (Q_e - Q_w)) ) ) - ...
             ( ((theta*b_dot_P(2,:)) + ((1-theta)*b_dot_P(1,:))) .* dx_P );
% in meters^2/yr
  
  
  residual(1) = ((h_P_t(1) - h_P(1)) .* ((dx_P(1))./dt_P)) + ...
                ( (1./W_P(1)) .* ( (theta* (Q_e_t(1) - Q_out_L_t)) + ...
                              ((1-theta)* (Q_e(1) - Q_out_L)) ) ) - ...
                ( ((theta*b_dot_P(2,1)) + ((1-theta)*b_dot_P(1,1))) .* (dx_P(1)) );

  
  residual(end) = ((h_P_t(end) - h_P(end)) .* ((dx_P(end))./dt_P)) + ...
                ( (1./W_P(end)) .* ( (theta* (Q_out_R_t - Q_w_t(end))) + ...
                              ((1-theta)* (Q_out_R - Q_w(end))) ) ) - ...
                ( ((theta*b_dot_P(2,end)) + ((1-theta)*b_dot_P(1,end))) .* (dx_P(end)) );
            
            
            
% -------------------------------------------------------------------------
if (full_domain_flag == 1)
% -------------------------------------------------------------------------


  Q_w(1)     = 0;    % by definition for the boundary condition of full domain
  Q_e(end)   = 0;     
  


% mass conservation in the wedges.
% ===============================
% calculate wedge length using new thickness values
% (m = slope; b = y-intercept; L = wedge length)
    m_wedge_R        = (h_P(end) - h_P(end-1)) / (x_P(end) - x_P(end-1));
    b_wedge_R        = (- m_wedge_R * x_P(end)) + h_P(end);
    wedge_length_R   = - ( b_wedge_R / m_wedge_R) - (x_w(end));

    m_wedge_R_t      = (h_P_t_guess(end) - h_P_t_guess(end-1)) / (x_P(end) - x_P(end-1));
    b_wedge_R_t      = (- m_wedge_R_t * x_P(end)) + h_P_t(end);
    wedge_length_R_t = - ( b_wedge_R_t / m_wedge_R_t) - (x_w(end));
    
  
    m_wedge_L        = (h_P(2) - h_P(1)) / (x_P(2) - x_P(1));
    b_wedge_L        = (- m_wedge_L * x_P(1)) + h_P(1);
    wedge_length_L   = abs(- ( b_wedge_L / m_wedge_L) - (x_e(1)));

    m_wedge_L_t      = (h_P_t_guess(2) - h_P_t_guess(1)) / (x_P(2) - x_P(1));
    b_wedge_L_t      = (- m_wedge_L_t * x_P(1)) + h_P_t(1);
    wedge_length_L_t = abs(- ( b_wedge_L_t / m_wedge_L_t) - (x_e(1)));
    
    
   


    h_w_use   = (h_P(end) + h_P(end-1))/2;
    h_w_t_use = (h_P_t(end) + h_P_t(end-1))/2;

    h_e_use   = (h_P(1) + h_P(2))/2;
    h_e_t_use = (h_P_t(1) + h_P_t(2))/2;
    
    
    Q_w_t(end) = -K_w_t(end) * ( (h_P_t(end) - h_P_t(end-1)) / ...
                                (x_P(end) - x_P(end-1)) );
                            
    Q_e_t(1)   = -K_e_t(1) * ( (h_P_t(2) - h_P_t(1)) / ...
                                (x_P(2) - x_P(1)) );
                            
                            
    residual(1)   = (0.5 * h_e_t_use * wedge_length_L_t * W_e(1) * (1/dt_P)) - ...
                    (0.5 * h_e_use * wedge_length_L * W_e(1) * (1/dt_P)) + ...
                    ( (theta * Q_e_t(1)) + ((1-theta) * Q_e(1)) ) - ...
                    ( (theta * b_dot_P(2,1) * wedge_length_L_t * W_e(1)) + ...
                      ((1-theta) * b_dot_P(1,1) * wedge_length_L * W_e(1)) );
    
                  
    residual(end) = (0.5 * h_w_t_use * wedge_length_R_t * W_w(end) * (1/dt_P)) - ...
                    (0.5 * h_w_use * wedge_length_R * W_w(end) * (1/dt_P)) - ...
                    ( (theta * Q_w_t(end)) + ((1-theta) * Q_w(end)) ) - ...
                    ( (theta * b_dot_P(2,end) * wedge_length_R_t * W_w(end)) + ...
                      ((1-theta) * b_dot_P(1,end) * wedge_length_R * W_w(end)) );
          
                  
    residual(1)   = residual(1) / W_e(1);                    
    residual(end) = residual(end) / W_w(end);              
            

end

          