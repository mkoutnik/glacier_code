function  [ Q_0_in, ...
            Q_ext_nodes, ...
            Q_ext ] = load_Q( t_nodes, t_P )


% ----------------------------------------------------------------------- %

%  get boundary flux at left edge of grid x_nodes(1) at times t_nodes
%  interpolate onto mesh.

% ----------------------------------------------------------------------- %


% Guess the starting value of flux entering at the left boundary
% of the limited domain
% ======================


  Q_0_in = (-10*800*(4/5));  % *25000);
%   Q_0_in = 0;

 
% Prescribe the value of external flux at left (L)
% and right (R) ends of limited domain
% =====================================

% ZERO External flux:
% -------------------
  Q_ext_L_nodes = 0 * ones(size(t_nodes));
  Q_ext_R_nodes = 0 * ones(size(t_nodes));

  
% % Linearly decreasing flux:
% % -------------------------
%   Q_ext_L_nodes = 0 * ones(size(t_nodes));
%   Q_ext_R_nodes = interp1([t_nodes(1) t_nodes(end)], [0 100000], t_nodes);
      
  

  Q_ext_nodes = [ Q_ext_L_nodes Q_ext_R_nodes ];   % N_t_nodes x 2
  

  
% interpolate Q_ext onto t_P
% ===========================
  Q_ext_L = interp1( t_nodes, Q_ext_L_nodes, t_P );
  Q_ext_R = interp1( t_nodes, Q_ext_R_nodes, t_P );
  

  Q_ext = [ Q_ext_L Q_ext_R ];   % N_t_nodes x 2
  
  
  
  
