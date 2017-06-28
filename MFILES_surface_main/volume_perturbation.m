function [ vol_pert ] = volume_perturbation ( x_flux, ...
                                              x_edges, ...
                                              dx_P, W, ...
                                              bdot_edges, ...
                                              bdot_edges_t )



% ------------------------------------------------------------------
% M. Koutnik, 29 Oct 2006

% adapted from flux_kin.m, written by Ed Waddington
% volume perturbation is convolved with the impulse response function
% the vol perturbation is the right hand side of:

% Q(x_flux(j), t_flux ) - Q(x_flux(j), t_flux-1) = 
%       Q_in( t_flux ) - Q_in( t_flux-1 )
%       + integral_{x_div(t)}^{x_flux(j)}[ b_dot(x,t) - b_dot(x,t-1)] W(x)dx
 
% This would force any additional mass out of the domain in one timestep
% Convolve this with the impulse response function to know how much
% mass should be exported so things are glaciologically correct!

% -------------------------------------------------------------------------

              
%  get intervals dx(i) = dx_P between x_w(i) and x_e(i) = x_w(i+1)
%  ===============================================================
       dx = dx_P;
       

%  add an extra interval to accommodate requests for flux at x(end)
%  ================================================================
       dx = [ dx dx(end) ];


       
%   ------------------------------------
%   SET UP NONDIMENSIONAL ARRAY X(i,j)
%   ------------------------------------

   %   X(i,j) = (x_ij - x_i) / (x_i+1 - x_1) = (x_ij- x_i) / dx_i 
   %   column j gives the  nondimensional value X of the point at which the
   %   integration ends within each x_node-defined interval [ x_i  x_i+1 ]
   %  each x_S interval   [ x_i  x_i+1 ] 
%  =======================================

      X = find_X_intervals( x_flux, x_edges );


% ------------------------
%  CREATE MATRIX bdot_t
% ------------------------

%  bdot_t contains the flux increment in each interval
%  by evaluating the analytical integral of bdot_edges * W over the interval
%  ========================================================================


% replace repmat calls here, this is called many times.
% A = s(ones(m,n));          % Equivalent to A = repmat(s,m,n);
% A = y(:,ones(1,n));        % Equivalent to A = repmat(y,1,n);
% A = x(ones(1,m),:);        % Equivalent to A = repmat(x,m,1)


   %   set up matrix b_i with bdot_edges(i) i=1:N_S  in each of length(x) columns
   %   Each column of b_i_plus1 is bdot_edges(i+1)  i=1:N_S
   %   ---------------------------------------------------
%          b_i_t       = repmat( bdot_edges_t', 1, length(x_flux) );   
           b_i_t       = bdot_edges_t( ones(1, length(x_flux)), : )';


         b_i_plus1_t = repmat( [bdot_edges_t(2:end) bdot_edges_t(end)]', ...
                                1, length(x_flux) );

         
   %  likewise matrix W_i columns are W(i)
   %   Each column of W_i_plus1 is W(i+1)
   %  ------------------------------------
%         W_i       = repmat( W', 1, length(x_flux) );
          W_i       = W ( ones(1, length(x_flux)), : )';
       
          
          W_i_plus1 = repmat( [ W(2:end) W(end)]', 1, length(x_flux) );

       

   % set up matrix with actual dx(i) intervals constituting each column
   %  -----------------------------------------------------------------
      %   dx_i  =  repmat( dx', 1, length(x_flux) );

        dx_i = dx( ones(1, length(x_flux)), :)';
      

  %  form integrals of products of 2 interpolation functions 
  %  for each interval i and for each x_flux(j)
  %  ==========================================
  
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
   %  ===================================================
   
         bdot_t = dx_i .* ( b_i_t .* W_i .* int_g0_g0 +  ...
                  ( b_i_t .*W_i_plus1 + b_i_plus1_t .* W_i )  .* int_g0_g1 ...
                  + b_i_plus1_t .* W_i_plus1 .* int_g1_g1 );



% ---------------------------
%  CREATE MATRIX bdot_tminus1
% ---------------------------

%  bdot_tminus1 contains the flux increment in each interval
%  by evaluating the analytical integral of bdot_edges * W over the interval
%  ========================================================================


   %   set up matrix b_i with bdot_edges(i) i=1:N_S  in each of length(x) columns
   %   Each column of b_i_plus1 is bdot_edges(i+1)  i=1:N_S
   %   ---------------------------------------------------
     %     b_i_tminus1       = repmat( bdot_edges', 1, length(x_flux) );
          b_i_tminus1 = bdot_edges( ones(1, length(x_flux)), :)';
         

       b_i_plus1_tminus1 = repmat( [bdot_edges(2:end) bdot_edges(end)]', 1, length(x_flux) );

         

         
%    % set up matrix with actual dx(i) intervals constituting each column
%    %  -----------------------------------------------------------------
%          dx_i  =  repmat( dx', 1, length(x_flux) );


  %  form integrals of products of 2 interpolation functions 
  %  for each interval i and for each x_flux(j)
  %  ==========================================
  
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
         bdot_tminus1 = dx_i .* ( b_i_tminus1 .* W_i .* int_g0_g0 +  ...
                       ( b_i_tminus1 .*W_i_plus1 + b_i_plus1_tminus1 .* W_i )  .* int_g0_g1 ...
                       + b_i_plus1_tminus1 .* W_i_plus1 .* int_g1_g1 );


% -----------
%  FIND flux
% -----------

  % no change in Qin here, because accounting for divide perturbation
  % later.
         vol_pert = (sum( bdot_t ) - sum( bdot_tminus1 ));

         
   

%  clean up the trash
%  ------------------
     clear  b_i_t b_i_tminus1  b_i_plus1_t b_i_plus1_tminus1  
     clear  W_i  W_i_plus1 bdot_t bdot_t_minus1

     