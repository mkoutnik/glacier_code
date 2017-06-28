function [ full_x_P, full_x_w, ...
           full_x_e, full_dx_P, ...
           full_dx_w, full_dx_e, ...
           full_h_P, full_h_w, full_h_e, ...
           full_S_P, full_S_w, full_S_e, ...
           wedge_length_L, ...
           wedge_length_R, ...
           ablation_rate_L, ...
           ablation_rate_R ] = wedge_terminus( full_x_P_in, full_x_w_in, ...
                                               full_x_e_in, full_dx_P_in, ...
                                               full_dx_w_in, full_dx_e_in, ...
                                               full_h_P_in, full_h_w_in, ...
                                               full_h_e_in, full_S_P_in, ...
                                               full_S_w_in, full_S_e_in, ...
                                               full_B_w_in, full_B_e_in, ...
                                               Q_in_wedge_L, ...
                                               Q_in_wedge_R, ...
                                               terminus_cutoff_L, ...
                                               terminus_cutoff_R )


% ----------------------------------------------------------------------- %
% Michelle Koutnik
% last updated: January 2009
%
% Detailed comments at the END of this file.
%
% ----------------------------------------------------------------------- %




% left boundary:
% ==============
  wedge_slope_L   = ( full_S_P_in(terminus_cutoff_L+1) - ...
                      full_S_P_in(terminus_cutoff_L) ) / ...
                      ( full_x_P_in(terminus_cutoff_L+1) - ...
                        full_x_P_in(terminus_cutoff_L) );

  wedge_slope_L   = - wedge_slope_L;   % switch it here for the calculation of L.
  
  
  y_intercept_L   = (- wedge_slope_L * abs(full_x_P_in(terminus_cutoff_L)) ) + ...
                     full_S_P_in(terminus_cutoff_L);

  total_length_L  = (full_B_e_in(terminus_cutoff_L) - y_intercept_L) / wedge_slope_L;
  wedge_length_L  = ( abs(total_length_L) - abs(full_x_e_in(terminus_cutoff_L)) );  % wedge starts at a western edge
  ablation_rate_L = - abs(Q_in_wedge_L / wedge_length_L);



% right boundary:
% ===============
  wedge_slope_R   = ( full_S_P_in(terminus_cutoff_R) - ...
                      full_S_P_in(terminus_cutoff_R-1) ) / ...
                      ( full_x_P_in(terminus_cutoff_R) - ...
                        full_x_P_in(terminus_cutoff_R-1) );

  y_intercept_R   = (- wedge_slope_R * full_x_P_in(terminus_cutoff_R)) + ...
                     full_S_P_in(terminus_cutoff_R);

  total_length_R  = (full_B_w_in(terminus_cutoff_R) - y_intercept_R) / wedge_slope_R;
  wedge_length_R  = abs( total_length_R - (full_x_w_in(terminus_cutoff_R)) );  % wedge starts at a western edge

  ablation_rate_R = - abs(Q_in_wedge_R / wedge_length_R);


% In resetting these values, don't forget:
%   - dx_w(1) includes part of the domain off the grid
%   - dx_e(end) includes part of the domain off the grid
%   - first boundary ends at a western edge, last boundary ends at an
%     eastern edge (!); both boundaries are cutoff at center points


% center points are just values from terminus cutoff on both sides
  full_x_P = full_x_P_in(terminus_cutoff_L: terminus_cutoff_R);
  full_h_P = full_h_P_in(terminus_cutoff_L: terminus_cutoff_R);  
  
  full_S_P = full_S_P_in(terminus_cutoff_L: terminus_cutoff_R);  
  
  
% western edge on left is zero thickness, 
% eastern edge on right is zero thickness
  full_x_w = [ sign(full_x_w_in(terminus_cutoff_L))*total_length_L ...
               full_x_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
         
  full_x_e = [ full_x_e_in(terminus_cutoff_L: terminus_cutoff_R-1) ...
               total_length_R ];

  full_h_w = [ 0 full_h_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
  full_h_e = [ full_h_e_in(terminus_cutoff_L: terminus_cutoff_R-1) 0 ];

  full_S_w = [ full_B_w_in(1) full_S_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
  full_S_e = [ full_S_e_in(terminus_cutoff_L: terminus_cutoff_R-1) full_B_e_in(end) ];
  
  
% now the differences...           
% dx_P is lengths of finite volumes
  full_dx_P = [ wedge_length_L ...
                full_dx_P_in(terminus_cutoff_L+1: terminus_cutoff_R-1) ...
                wedge_length_R ];
          
% first dx_w is off the grid, so only keep part of it.            
  full_dx_w = [ wedge_length_L - (full_dx_w_in(terminus_cutoff_L)/2) ...
                full_dx_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
                
  full_dx_e = [ full_dx_e_in(terminus_cutoff_L: terminus_cutoff_R-1) ...
               wedge_length_R - (full_dx_e_in(terminus_cutoff_R-1)/2) ];

           


% CHECK WITH A PLOT:
% ==================

% plot(full_x_P_in, full_S_P_in)
% hold on
% plot(full_x_P_in, full_S_P_in,'b.')

%plot([full_x_w_in(1) full_x_e_in], [full_B_w_in(1) full_B_e_in], 'r')
%plot( [-total_length_L 0], [full_B_w_in(1) y_intercept_L],'c')    
%plot( [total_length_R 0], [full_B_e_in(end) y_intercept_R],'c')    













% % left boundary:
% % ==============
%   wedge_slope_L   = ( full_S_P_in(terminus_cutoff_L+1) - ...
%                       full_S_P_in(terminus_cutoff_L) ) / ...
%                       ( full_x_P_in(terminus_cutoff_L+1) - ...
%                         full_x_P_in(terminus_cutoff_L) );
% 
%   wedge_slope_L   = - wedge_slope_L;   % switch it here for the calculation of L.
%   
%   
%   y_intercept_L   = (- wedge_slope_L * abs(full_x_P_in(terminus_cutoff_L)) ) + ...
%                      full_S_P_in(terminus_cutoff_L);
% 
%   total_length_L  = (full_B_e_in(terminus_cutoff_L) - y_intercept_L) / wedge_slope_L;
%   wedge_length_L  = ( abs(total_length_L) - abs(full_x_e_in(terminus_cutoff_L)) );  % wedge starts at a western edge
%   
%   wedge_length_L_for_area = sqrt(wedge_length_L^2 + (full_S_e_in(terminus_cutoff_L)-full_B_e_in(terminus_cutoff_L))^2);
%   
%   ablation_rate_L = - abs(Q_in_wedge_L / (wedge_length_L_for_area));
% 
%   % ablation_rate_L = - abs(Q_in_wedge_L / (wedge_length_L * (full_S_P_in(terminus_cutoff_L)-full_B_e_in(terminus_cutoff_L))));
% 
% 
% 
% % right boundary:
% % ===============
%   wedge_slope_R   = ( full_S_P_in(terminus_cutoff_R) - ...
%                       full_S_P_in(terminus_cutoff_R-1) ) / ...
%                       ( full_x_P_in(terminus_cutoff_R) - ...
%                         full_x_P_in(terminus_cutoff_R-1) );
% 
%   y_intercept_R   = (- wedge_slope_R * full_x_P_in(terminus_cutoff_R)) + ...
%                      full_S_P_in(terminus_cutoff_R);
% 
%   total_length_R  = (full_B_w_in(terminus_cutoff_R) - y_intercept_R) / wedge_slope_R;
%   wedge_length_R  = abs( total_length_R - (full_x_w_in(terminus_cutoff_R)) );  % wedge starts at a western edge
% 
%   wedge_length_R_for_area = sqrt(wedge_length_R^2 + (full_S_w_in(terminus_cutoff_R)-full_B_w_in(terminus_cutoff_R))^2);
%   
%   ablation_rate_R = - abs(Q_in_wedge_R / wedge_length_R_for_area);
% 
% 
% % In resetting these values, don't forget:
% %   - dx_w(1) includes part of the domain off the grid
% %   - dx_e(end) includes part of the domain off the grid
% %   - first boundary ends at a western edge, last boundary ends at an
% %     eastern edge (!); both boundaries are cutoff at center points
% 
% 
% % center points are just values from terminus cutoff on both sides
%   full_x_P = full_x_P_in(terminus_cutoff_L: terminus_cutoff_R);
%   full_h_P = full_h_P_in(terminus_cutoff_L: terminus_cutoff_R);  
%   
%   full_S_P = full_S_P_in(terminus_cutoff_L: terminus_cutoff_R);  
%   
%   
% % western edge on left is zero thickness, 
% % eastern edge on right is zero thickness
%   full_x_w = [ sign(full_x_w_in(terminus_cutoff_L))*total_length_L ...
%                full_x_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
%          
%   full_x_e = [ full_x_e_in(terminus_cutoff_L: terminus_cutoff_R-1) ...
%                total_length_R ];
% 
%   full_h_w = [ 0 full_h_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
%   full_h_e = [ full_h_e_in(terminus_cutoff_L: terminus_cutoff_R-1) 0 ];
% 
%   full_S_w = [ full_B_w_in(1) full_S_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
%   full_S_e = [ full_S_e_in(terminus_cutoff_L: terminus_cutoff_R-1) full_B_e_in(end) ];
%   
%   
% % now the differences...           
% % dx_P is lengths of finite volumes
%   full_dx_P = [ wedge_length_L ...
%                 full_dx_P_in(terminus_cutoff_L+1: terminus_cutoff_R-1) ...
%                 wedge_length_R ];
%           
% % first dx_w is off the grid, so only keep part of it.            
%   full_dx_w = [ wedge_length_L - (full_dx_w_in(terminus_cutoff_L)/2) ...
%                 full_dx_w_in(terminus_cutoff_L+1: terminus_cutoff_R) ];
%                 
%   full_dx_e = [ full_dx_e_in(terminus_cutoff_L: terminus_cutoff_R-1) ...
%                wedge_length_R - (full_dx_e_in(terminus_cutoff_R-1)/2) ];
% 
%            
% 
% 
% % % CHECK WITH A PLOT:
% % % ==================
% 
% % plot(full_x_P_in, full_S_P_in)
% % plot(full_x_P_in, full_S_P_in,'b.')
% 
% % plot([full_x_w_in(1) full_x_e_in], [full_B_w_in(1) full_B_e_in], 'r')
% % plot( [-total_length_L 0], [full_B_w_in(1) y_intercept_L],'c')    
% % plot( [total_length_R 0], [full_B_e_in(end) y_intercept_R],'c')    





% 
% % %plot(full_x_P_in, full_S_P_in)
% % plot(full_x_P_in, full_h_P_in)
% % hold on
% % %plot(full_x_P_in, full_S_P_in,'b.')
% % plot(full_x_P_in, full_h_P_in,'b.')
% % 
% % plot(full_x_P_in, zeros(size(full_x_P_in)), 'r')
% % plot( [-total_length_L 0], [0 y_intercept_L],'c')    
% % plot( [total_length_R 0], [0 y_intercept_R],'c')    
% % 
% % % plot([full_x_w_in(1) full_x_e_in], [full_B_w_in(1) full_B_e_in], 'r')
% % % plot( [-total_length_L 0], [full_B_w_in(1) y_intercept_L],'c')    
% % % plot( [total_length_R 0], [full_B_e_in(end) y_intercept_R],'c')    
% 
