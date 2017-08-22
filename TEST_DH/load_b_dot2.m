function [ b_dot_nodes, ...
           b_dot_P, ...
           b_dot_edges ] = load_b_dot( x_P, x_edges, t_P, ...
                                       x_nodes, t_nodes )

% -------------------------------------------------------------------------                                   
                                   
global N_t_nodes2 N_x_nodes2
  

global DIRECTORY_data

addpath(DIRECTORY_data)


% Interpolate onto nodes and mesh:
% ================================
  b_dot_nodes = NaN * ones(N_t_nodes2, N_x_nodes2);
  
  
  
% Estimate from available observations
% ====================================

% From QGIS -- Accumulation_A vs. Accumulation_R
load DH_accum_width_velocity.mat

load DH_surf_bed.mat

% if you want to calculate from a lapse rate similar to Bliss et al. (2011) for Taylor glacier 
precip_at_sl = -0.35;
lapse = 0.35/1500;

Hat_bdot_modern_lapse = precip_at_sl + lapse.*Hat_modern_surface;

b_dot_use = interp1(Hat_centerline_distance, Hat_bdot_modern_lapse, x_nodes);

% % b_dot_use = interp1(Hat_accumulation_A_centerline_distance, Hat_accumulation_A, x_nodes);
% b_dot_use = -interp1(Hat_accumulation_centerline_distance, Hat_accumulation_A, x_nodes);
 

% % Use Trevor's estimate: I just calculated it from a lapse rate and the modern surface DEM:
% % precip_at_sl = -0.35;  %precipitation at sea level
% % lapse = 0.34/1500; %precipitation lapse rate with elevation (m/m)
% % Darwin_bdot_modern= precip_at_sl + lapse.*Darwin_modern_surface
% % I took these values from Bliss et al. (2011) at Taylor Glacier, and tweaked them slightly to match the surface of Darwin. (0.34/1500 instead of 0.35/1500)
% 
% load Darwin_bdot_modern.mat
% %b_dot_use = interp1(Darwin_x_P, Darwin_bdot_modern, x_nodes, 'linear', 'extrap');
% b_dot_use = interp1([Darwin_x_P x_nodes(end)], [Darwin_bdot_modern Darwin_bdot_modern(end)], x_nodes);




  for ii = 1:N_t_nodes2
    b_dot_nodes(ii,:) = b_dot_use;
  end

  

% % Constant in time
% % ------------------
%  b_dot_use = interp1([x_P(1) x_P(end)], [0.3 0.3], x_P);  % Average accumulation 28 cm/yr
%                                                           % Denton et al. (1989)                                                                                                          
%  for ii = 1:N_t_nodes
%    b_dot_nodes(ii,:) = interp1(x_P, b_dot_use, ...
%                                x_nodes, 'linear', 'extrap');
%  end

       
% % Linearly decreasing in time
% % ----------------------------
% b_dot_use      = interp1([x_P(1) x_P(end)], [0.1 0.1], x_P);
% time_variation = interp1([t_P(1) t_P(end)], [6 3], t_nodes); 
% 
% for ii = 1:N_t_nodes
%     b_dot_nodes(ii,:) = time_variation(ii) * interp1(x_P, b_dot_use, x_nodes, 'linear', 'extrap');
% end
 




% interpolate onto mesh -- from nodes!
% ----------------------                                        
   b_dot_P = interp2( repmat( x_nodes, N_t_nodes2, 1 ), ...
                      repmat( t_nodes, 1, N_x_nodes2 ), ...
                      b_dot_nodes, x_P, t_P );

 
% interpolate onto edges of mesh
% ------------------------------ 
   b_dot_edges = interp2( repmat( x_nodes, N_t_nodes2, 1 ), ...
                           repmat( t_nodes, 1, N_x_nodes2 ), ...
                           b_dot_nodes, x_edges, t_P );                                     
                                               