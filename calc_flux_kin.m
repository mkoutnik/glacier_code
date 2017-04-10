   function flux = calc_flux_kin( x_flux, x_edges, dx_P, W, ...
                                  S_dot, bdot_edges, Q_in_t_flux )


% -------------------------------------------------------------------------
% from Ed Waddington

% calculate kinematic flux.
% Detailed comments at the end of this file.
% -------------------------------------------------------------------------

                       

% %  set up widths W at all control-volumes edges
%        W = [ W_w(1)  W_e ];
% % 
% 
% %  Q_in_t_flux  is current boundary flux
% %  S_dot_t_flux is current rate of thickness change
% %  bdot_edges(x)  contains current accumulation rates at 
% %               all control-volume edges
% %  -------------------------------------------------------------
%        bdot_edges  = [ b_dot_w(1)  b_dot_e ];
% %

%  x_edges contains locations of all control-volumes edges
%       x_edges  = [ x_w(1) x_e ];
%

%  get intervals dx(i) = dx_P between x_w(i) and x_e(i) = x_w(i+1)
%  ---------------------------------------------------------------
       dx_temp = dx_P;
%
   %  add an extra interval to accommodate requests for flux at x(end)
       dx = [ dx_temp ]; % dx_temp(end) ];


%   ------------------------------------
%   SET UP NONDIMENSIONAL ARRAY X(i,j)
%   ------------------------------------
   %   X(i,j) = (x_ij - x_i) / (x_i+1 - x_1) = (x_ij- x_i) / dx_i 
   %   column j gives the  nondimensional value X of the point at which the
   %   integration ends within each x_node-defined interval [ x_i  x_i+1 ]
   %  each x_S interval   [ x_i  x_i+1 ] 
%  ---------------------------------------------------------------

      X = find_X_intervals( x_flux, x_edges );


% ------------------------
%  CREATE MATRIX Q_bits
% ------------------------
%  Q_bits contains the flux increment in each interval
%  by evaluating the analytical integral of bdot_edges * W over the interval
%  --------------------------------------------------------------------

   %   set up matrix b_i with bdot_edges(i) i=1:N_S  in each of length(x) columns
   %   Each column of b_i_plus1 is bdot_edges(i+1)  i=1:N_S
   %   ---------------------------------------------------
         b_i       = repmat( bdot_edges', 1, length(x_flux) );
         b_i_plus1 = repmat( [bdot_edges(2:end) bdot_edges(end)]', 1, length(x_flux) );
 
   %  likewise matrix W_i columns are W(i)
   %   Each column of W_i_plus1 is W(i+1)
   %  ------------------------------------
        W_i       = repmat( W', 1, length(x_flux) );
        W_i_plus1 = repmat( [ W(2:end) W(end)]', 1, length(x_flux) );


   % set up matrix with actual dx(i) intervals constituting each column
   %  -----------------------------------------------------------------
         dx_i  =  repmat( dx', 1, length(x_flux) );


  %  form integrals of products of 2 interpolation functions 
  %  for each interval i and for each x_flux(j)
  %  ------------------------------------------
       %   gamma_0^2(X)     -  gamma_0 = X   on [  0 <= X <= 1 ]
       %   -----------------------------------------------------
             int_g0_g0 = ( 1 - ( 1 - X).^3  ) /3;

       %   gamma_1^2(X)     -  gamma_1 = ( 1 - X )   on [  0 <= X <= 1 ]
       %   -------------------------------------------------
             int_g1_g1 = X.^3 /3;

       %   gamma_0(X) .* gamma_1(X)   
       %   ---------------------------
             int_g0_g1 = X.^2 .* ( 1/2 - X/3);


   %  Q_bits(i,j) is the integral of bdot_edges * W in x_S interval "i"
   %  as appropriate for finding flux at x_flux(j)
   %  (intervals beyond x_flux(j) make zero contribution)
   %  ---------------------------------------------------
         Q_bits = dx_i .* ( b_i .* W_i .* int_g0_g0 +  ...
                  ( b_i .*W_i_plus1 + b_i_plus1 .* W_i )  .* int_g0_g1 ...
                  + b_i_plus1 .* W_i_plus1 .* int_g1_g1 );

%
%
% --------------------------
%  CREATE MATRIX Store_bits
% --------------------------

%  Store_bits contains the flux increment in each interval
%  by evaluating the analytical integral of S_dot * W over the interval
%  --------------------------------------------------------------------

   % integrals of basis functions between zero and X
      int_0 = X.*(1-X/2);
      int_1 = X.*X/2;

   % stuff matrix with a row of S_dot for each x_flux(i)
     % S_dot_i = repmat( [ S_dot S_dot(end) ]', 1, length(x_flux) );
      S_dot_i = repmat( [ S_dot ]', 1, length(x_flux) );
      
%  multiply to get storage in each control volume, for each x_flux(i)
      store_bits = dx_i .* S_dot_i .*(  (W_i .* int_0) + (W_i_plus1 .* int_1) );


% -------------
%  FIND flux
% -------------
   %  sum columns and add Q_in_t_flux to find total flux at x_j
   %  ---------------------------------------------------------
         flux = Q_in_t_flux + sum( Q_bits ) - sum( store_bits );
  
         
%  clean up the trash
%  ------------------
 %    clear  b_i  b_i_plus1  W_i  W_i_plus1   Q_bits


     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate kinematic flux Q( x_flux, t_flux )  -  row vector
%  x_flux must be a row vector of positions
%        Here, t_flux is a scalar
%
%     Q(x_flux(j), t_flux ) = Q_in( t_flux )
%                             + integral_{x_eiv(t)}^{x_flux(j)} 
%                                [ b_dot(x,t) - S_dot(x,t)] W(x)  dx 
%
%   x_i are the control-volume edges x_edge = [ x_w(1)  x_eS ] at which
%     b_dot and W are specified
%   x_flux(j) are the points at which we want to find Q( x_flux(j), t_flux )
%     at time t_flux
%   x_ij are the points at which integration ends in interval i =  [ x_i  x_i+1 ]
%   when finding  Q(x_j).    
%   x_ij = x_i+1  when x_i+1 < x_flux_j
%             = x_j     when x_i < x_flux_j  < x_i+1
%             = x_i     when x_flux_j < x_i
%
%  b_dot(x,t) and W(x) are linear between x_edge(i) and x_edge(i+1)
%        
%   First get vector bdot_edges(x) by interpolating b_dot(x,t) at time t_flux  
%      bdot_edges(i) contains accumulation rates at x_edge(i) at time t_flux
%
%        bdot_edges(x) = bdot_edges_i gamma_{0i}(X) + bdot_edges{i+1} gamma_{1i}(X)
%        W(x)        = W_i gamma_{0i}(X)     + W_{i+1} gamma_{1i}(X)
%
%   where X = (x - x_i)/(x_i+1 - x_i)
%   and linear interpolating functions are
%      gamma_{0i}(X) = 1 at X=0
%                    = 0 at X=1
%      gamma_{1i}(X) = 0 at X=0
%                    = 1 at X=1
%
%  All intervals except the last span a control-volume,  
%  and as a result have a simple integration formula.
%  Last interval ends at an arbitrary point inside a control volume,
%  so its formula is more complicated.
%  The general formula for integral between x_i and x_i+1
%
%  dQ = (dx_i/3) ( bdot_edges_i {  W_i ((x-x_i)/dx_i)^3                          }
%                          {+ W_i+1 ((x-x_i)/dx_i)^2 (3/2- ((x-x_i)/dx_i) ) }
%             + bdot_edges_i+1 { W_i+1 (1 - ((x_i+1 - x)/dx_i)^3 )            }
%                         {+ W_i ((x-x_i)/dx_i)^2 (3/2- ((x-x_i)/dx_i) ) } ) 
%
%     x_i are the control-volume edges where W(x) and b_dot(x) are defined
%     and dx_i = (x_i+1 - x_i)
%
% -------------------------------------------------------------------------


%
%

