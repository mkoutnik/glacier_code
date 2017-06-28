function [ x_ldm_P, h_ldm_P, x_ldm_div, ...
           x_div_full, index_div_full, ...
           full_x_P, full_x_w, full_x_e, ...
           full_dx_P, full_dx_w, full_dx_e, ...
           full_h_P, full_h_w, full_h_e, ...
           full_S_P, full_S_w, full_S_e, ...
           full_B_P, full_B_w, full_B_e, ...
           full_W_P, full_W_w, full_W_e, ...
           full_slip_P, full_slip_w, full_slip_e, ...
           full_b_dot_P, full_b_dot_w, full_b_dot_e, ...
           full_A_eff_edges, ...
           wedge_length_L, ...
           wedge_length_R, ...
           Q_in_wedge_L, ...
           Q_in_wedge_R, ...
           ablation_rate_L, ...
           ablation_rate_R ] = ldm_extension ( h_ldm_P, x_P, ...
                                               x_e, x_w, ...
                                               dx_P, dx_e, dx_w, ...
                                               B_P, W_P, slip_P, ...
                                               A_eff_edges, ...
                                               b_dot_P, x_div )
                       
                                           
% ----------------------------------------------------------------------- %
% Michelle Koutnik
% last updated: March 2010, Copenhagen


% Embed the limited model in a full-domain model:
%   1) estimate the length using a Paterson model
%   2) regrid for a nonuniform grid
%   3) setup a simple mass-balance pattern over the full length
%   4) using numerical surface solver to calculate full surface
%   5) replace terminus with a wedge shape



% LDM    = limited-domain model.
% '_ext' = represents values for the extended domain

%
% ----------------------------------------------------------------------- %


  global paterson_c_over_a
  global n
 
  
 

%   % try resetting bdot for the case when there is a strong gradient...
%   b_dot_P_orig = b_dot_P;
%   b_dot_P_mean = mean(b_dot_P_orig(1,:));
%   
%   b_dot_P = ones(size(b_dot_P)) * b_dot_P_mean;
  
  
%    A_eff_edges_orig = A_eff_edges;
%    A_eff_edges_mean = mean(A_eff_edges_orig);
%    
%    A_eff_edges = ones(size(A_eff_edges_orig)) * A_eff_edges_mean;




% -------------------------------------------------------------------------
% 1. Find full span length, L
% -------------------------------------------------------------------------
 
 
 disp( ' ' )
 disp(' Spacing for limited-domain calculation is reset for WAIS-size ice sheet! ')
 disp( ' ' )
 
 
 dx_ext_P   = mean(dx_P);
 coarse_dx  = 5000;
 value_cut  = 2;
 add_length = 50000;    % add this to length estimate because accumulation
                        % pattern may not be uniform
 
 
% % if using a thinner ice sheet
% dx_ext_P  = mean(dx_P);
% value_cut = 2;
% coarse_dx  = 300;
% add_length = 20000;
 
                        

% x_ldm  = x_P-position at end of LDM
% x_div  = x_P-position of the divide in LDM
% H_0    = ice-sheet thickness at the divide
% h_0    = ice-sheet thickness at x_ldm
% L_ext  = length of extended domain (from the divide)

% these values are all on x_P, for the limited domain:
  x_ldm_L       = x_P(1);      
  x_ldm_R       = x_P(end); 
  x_ldm_P       = x_P;
  x_ldm_div     = x_div;
  index_ldm_div = find(x_P == x_ldm_div);
  
  S_0   = h_ldm_P(index_ldm_div) + B_P(index_ldm_div);
  S_0_L = h_ldm_P(1) + B_P(1);
  S_0_R = h_ldm_P(end) + B_P(end);
  
% use these to estimate length from Paterson model  
% this will overestimate the length when Bed /= 0 m, but cutoff the
% imaginary values of the surface profile later.
   H_0_flat   = S_0;
   h_0_L_flat = S_0_L;
   h_0_R_flat = S_0_R;
   
  
% Paterson-model length
% ======================
  term1   = 1 - ((h_0_L_flat/H_0_flat)^(2+(2/n)));
  term2   = (1 + paterson_c_over_a)^(1/n);
  L_ext_L = abs(x_ldm_L - x_div) / ...
            ( (term1 / term2).^(1 / (1+(1/n))) );
  L_ext_L = - L_ext_L;
  
        
  term1   = 1 - ((h_0_R_flat/H_0_flat)^(2+(2/n)));
  term2   = (1 + paterson_c_over_a)^(1/n);
  L_ext_R = (x_ldm_R - x_div) / ...
            ( (term1 / term2).^(1 / (1+(1/n))) );   
   
   

% Extend the domain
% =================
  full_x_P_L    = [ 0: -dx_ext_P: L_ext_L ];    
  full_x_P_R    = [ 0: dx_ext_P: L_ext_R ];
  full_x_P      = [ fliplr(full_x_P_L(2:end)) full_x_P_R ];
 
  x_div_full     = 0;
  index_div_full = find(full_x_P == 0);          % divide defined to be at zero.
 
  dx_P_temp1    = diff( full_x_P );
  full_dx_w     = [ dx_P_temp1(1) dx_P_temp1 ];   % first value is outside domain
  full_dx_e     = [ dx_P_temp1 dx_P_temp1(end) ]; % last value is outside domain

% change it so that first and last edges are full length.  
  full_x_w      = full_x_P - full_dx_w/2;
 % full_x_w(1)   = L_ext_L;
  full_x_e      = full_x_P + full_dx_e/2;
 % full_x_e(end) = L_ext_R;
  
  full_dx_P     = full_x_e - full_x_w;

  full_x_edges  = [ full_x_w(1) full_x_e ];
  
  
 % original grid across the limited domain:
   index_ldm_L = find((full_x_P <= (x_P(1) - x_div)), 1, 'last');
   index_ldm_R = find((full_x_P <= (x_P(end) - x_div)), 1, 'last'); 

  
  
                          
   
  
% -------------------------------------------------------------------------
% 3. Setup the mass-balance pattern
% -------------------------------------------------------------------------
  

% Paterson model assumes zone of uniform accumulation and zone of uniform
% ablation separated by discontinuity at equilibrium line R

% sorry for "R" and "R" -- underscrore R ("_R") is for right side of the
% limited domain; "_L" is for left side of limited domain (across divide)

  c = paterson_c_over_a;     % ratio of accumulation rate to ablation rate
  a = 1;

% equilibrium-line position on both sides of the divide
  R_L_actual = L_ext_L * (a/(a+c) );
  R_R_actual = L_ext_R * (a/(a+c) );
  
  % adjust equilibrium line so that it is exactly on a grid boundary.
  R_L = full_x_P(find(full_x_P < R_L_actual,1,'last'));
  R_R = full_x_P(find(full_x_P < R_R_actual,1,'last'));
 
   
  full_x_P_L = full_x_P(1:index_div_full);
  full_x_P_R = full_x_P(index_div_full:end);
      
  
% %  find points in accumulation and ablation area
%    index_c_L =  full_x_P_L >= R_L ;  % change signs! x negative on the left.
%    index_a_L =  full_x_P_L < R_L ;
%    index_c_R =  full_x_P_R <= R_R ;
%    index_a_R =  full_x_P_R > R_R ;


index_c_L = full_x_P_L >= 0.85*L_ext_L;
index_a_L = full_x_P_L < 0.85*L_ext_L;
index_c_R = full_x_P_R <= 0.85*L_ext_R;
index_a_R = full_x_P_R > 0.85*L_ext_R;

  
%  set profile
     c_c = (1 + c/a).^(1/n);
     c_a = (1 + a/c).^(1/(n+1) );
     h_ext_L(index_c_L) = H_0_flat * ( 1 - c_c*((full_x_P_L(index_c_L))/L_ext_L).^(1+1/n) ).^(n/(2*n+2) );
     h_ext_L(index_a_L) = H_0_flat * ( c_a*(1 - (full_x_P_L(index_a_L))/L_ext_L) ).^(1/2);
    
     h_ext_R(index_c_R) = H_0_flat * ( 1 - c_c*((full_x_P_R(index_c_R))/L_ext_R).^(1+1/n) ).^(n/(2*n+2) );
     h_ext_R(index_a_R) = H_0_flat * ( c_a*(1 - (full_x_P_R(index_a_R))/L_ext_R) ).^(1/2);
    
     
     full_h_P_pat = [ h_ext_L(1:end-1) h_ext_R ];
   
     
  
%    plot(full_x_P, full_h_P_pat, 'k')
%    hold on
    
    
    
% accumulation rates from values at limited-domain edges     
  acc_rate_L   = b_dot_P(1,1);
  acc_rate_R   = b_dot_P(1,end);
  abl_rate_L   = acc_rate_L / paterson_c_over_a;
  abl_rate_R   = acc_rate_R / paterson_c_over_a;
  
  
  full_b_dot_L = [ ones(size(full_x_P_L(index_a_L))) * -abl_rate_L ...
                    ones(size(full_x_P_L(index_c_L))) * acc_rate_L ];
  full_b_dot_R = [ ones(size(full_x_P_R(index_c_R))) * acc_rate_R ...
                    ones(size(full_x_P_R(index_a_R))) * -abl_rate_R ];

   

  full_b_dot_P = [ full_b_dot_L(1:end-1) full_b_dot_R ];

  
  % replace with actual values across limited domain:
    full_b_dot_P(index_ldm_L:index_ldm_R) = b_dot_P(1,:);
  

  [ full_b_dot_w, ...
    full_b_dot_e ] = get_edge_values_quadratic ( full_b_dot_P, full_x_P, ...
                                                 full_x_w, full_x_e, ...
                                                 full_dx_P, full_dx_w, ...
                                                 full_dx_e );        
 
 % figure(2)
 %   plot([full_x_w(1) full_x_e], [full_b_dot_w(1) full_b_dot_e])
                                      

 % extend other values  
 % ==================== 

  extend_length_L = length(full_x_P(1:index_div_full)) - ...
                    length(x_P(1:index_ldm_div));
  extend_length_R = length(full_x_P(index_div_full:end)) - ...
                    length(x_P(index_ldm_div:end));

% bed:                 
  full_B_P = [ ones(1,extend_length_L)*B_P(1) B_P ones(1,extend_length_R)*B_P(end) ];
[ full_B_w, ...
  full_B_e ] = get_edge_values_quadratic ( full_B_P, full_x_P, ...
                                           full_x_w, full_x_e, ...
                                           full_dx_P, full_dx_w, full_dx_e );     
               
% width:           
  full_W_P = [ ones(1,extend_length_L)*W_P(1) W_P ones(1,extend_length_R)*W_P(end) ]; 
[ full_W_w, ...
  full_W_e ] = get_edge_values_quadratic ( full_W_P, full_x_P, ...
                                           full_x_w, full_x_e, ...
                                           full_dx_P, full_dx_w, full_dx_e );     
  
% slip: 
  full_slip_P = [ ones(1,extend_length_L)*slip_P(1) slip_P ones(1,extend_length_R)*slip_P(end) ];           
[ full_slip_w, ...
  full_slip_e ] = get_edge_values_quadratic ( full_slip_P, full_x_P, ...
                                              full_x_w, full_x_e, ...
                                              full_dx_P, full_dx_w, full_dx_e );     
  
  
% effective softness:                         
  full_A_eff_edges = [ ones(1,extend_length_L)* A_eff_edges(1) ...
                       A_eff_edges ...
                       ones(1,extend_length_R)* A_eff_edges(end) ];                   
              

 
                   
% -------------------------------------------------------------------------
% 3. calculate the full-domain surface
% -------------------------------------------------------------------------
 
 % use numerical surface solver "calc_h_0.m"

 
  Q_in_use = 0.001;     % start at the divide; can't use exactly = 0 
  S_in_use = S_0;       % elevation from limited model
  
  
  
% left side of the divide:
% ------------------------
x_P_left = abs(fliplr(full_x_P(1:index_div_full)));
x_w_left = [full_x_w(index_div_full) ...             % switch order x_e and x_w
                              abs(fliplr(full_x_e(1:index_div_full-1)))];
x_e_left = abs(fliplr(full_x_w(1:index_div_full)));
b_dot_edges_left =  fliplr([full_b_dot_w(1) ...
                             full_b_dot_e(1:index_div_full)]);                    
                          
A_eff_w_left = fliplr(full_A_eff_edges(1:index_div_full));
A_eff_e_left = fliplr(full_A_eff_edges(2:index_div_full+1));
  

  [ full_S_P_L, ...
    full_h_P_L ] = calc_h_0( x_P_left, x_w_left, x_e_left, ...
                             fliplr(full_dx_P(1:index_div_full)), ...
                             fliplr(full_B_P(1:index_div_full)), ...
                             fliplr(full_B_w(1:index_div_full)), ...  % e values first???
                             fliplr(full_B_e(1:index_div_full)), ...
                             fliplr(full_W_P(1:index_div_full)), ...
                             fliplr(full_W_w(1:index_div_full)), ...
                             fliplr(full_W_e(1:index_div_full)), ...
                             b_dot_edges_left, ...
                             fliplr(full_slip_w(1:index_div_full)), ...
                             fliplr(full_slip_e(1:index_div_full)), ...
                             Q_in_use, S_in_use, ...
                             A_eff_w_left, A_eff_e_left );
                         
                         
% overestimated the length, so some values will be imaginary
  cutoff_L = (find(full_h_P_L < 0, 1, 'first'))-value_cut;        
  
  if ( (isempty(cutoff_L) == 1) && ( isempty (find(diff(full_S_P_L)>0, 1)) == 0) )
      cutoff_L = find(diff(full_S_P_L)>0, 1);
  
  elseif (isempty(cutoff_L) == 1)
      cutoff_L = length(full_h_P_L);  % cutoff at least one.
  end
  
% make sure all the imaginary values were removed.  
  if (isreal(full_h_P_L(1:cutoff_L)) ~= 1)
      disp(' check full-domain extension in ldm_extension.m')
      stop;
  end
  
  
                                             
% right side of the divide:
% ------------------------
  [ full_S_P_R, ...
    full_h_P_R ] = calc_h_0( full_x_P(index_div_full:end), ...
                             full_x_w(index_div_full:end), ...
                             full_x_e(index_div_full:end), ...
                             full_dx_P(index_div_full:end), ...
                             full_B_P(index_div_full:end), ...
                             full_B_w(index_div_full:end), ...
                             full_B_e(index_div_full:end), ...
                             full_W_P(index_div_full:end), ...
                             full_W_w(index_div_full:end), ...
                             full_W_e(index_div_full:end), ...
                             [full_b_dot_w(index_div_full) ...
                              full_b_dot_e(index_div_full:end)], ...
                             full_slip_w(index_div_full:end), ...
                             full_slip_e(index_div_full:end), ...
                             Q_in_use, S_in_use, ...
                             full_A_eff_edges(index_div_full:end-1), ...
                             full_A_eff_edges(index_div_full+1:end) );                                            
                                                                                 
                                             
  cutoff_R = (find(full_h_P_R < 0, 1, 'first'))-value_cut;   
  
  if ( isreal(full_h_P_R(1:cutoff_R)) ~= 1 )   % check if any imaginary
      cutoff_R = find(imag(full_h_P_R) == 0, 1, 'last');
  end
  
  if ( (isempty(cutoff_R) == 1) && ( isempty (find(diff(full_S_P_R)>0, 1)) == 0) )
      cutoff_R = find(diff(full_S_P_R)>0, 1);
  
  elseif (isempty(cutoff_R) == 1)
      cutoff_R = length(full_h_P_R);   % cutoff at least one.
  end                                           
    
  if (isreal(full_h_P_R(1:cutoff_R)) ~= 1)
      disp(' check full-domain extension in ldm_extension.m')
      stop;
  end
  
  
% put together the full surface:  
  full_h_P = [ fliplr(full_h_P_L(1:cutoff_L)) full_h_P_R(2:cutoff_R) ];
  full_S_P = [ fliplr(full_S_P_L(1:cutoff_L)) full_S_P_R(2:cutoff_R) ];
  
    
% how much of domain was used:
  diff_L   = length(full_S_P_L) - cutoff_L;
  diff_R   = length(full_S_P_R) - cutoff_R;
  
  
% resize:
  full_x_P         = full_x_P(diff_L+1:end-diff_R);
  full_x_w         = full_x_w(1+diff_L:end-diff_R);
  full_x_e         = full_x_e(1+diff_L:end-diff_R);
  full_x_edges     = full_x_edges(1+diff_L:end-diff_R);
  full_dx_P        = full_dx_P(1+diff_L:end-diff_R);
  full_dx_w        = full_dx_w(1+diff_L:end-diff_R);
  full_dx_e        = full_dx_e(1+diff_L:end-diff_R);
  full_W_P         = full_W_P(1+diff_L:end-diff_R);
  full_W_w         = full_W_w(1+diff_L:end-diff_R);
  full_W_e         = full_W_e(1+diff_L:end-diff_R);
  full_slip_P      = full_slip_P(1+diff_L:end-diff_R);
  full_slip_w      = full_slip_w(1+diff_L:end-diff_R);
  full_slip_e      = full_slip_e(1+diff_L:end-diff_R);
  full_A_eff_edges = full_A_eff_edges(1+diff_L:end-diff_R);
  full_b_dot_P     = full_b_dot_P(1+diff_L:end-diff_R);
  full_b_dot_w     = full_b_dot_w(1+diff_L:end-diff_R);
  full_b_dot_e     = full_b_dot_e(1+diff_L:end-diff_R);
  full_B_P         = full_B_P(1+diff_L:end-diff_R);
  full_B_w         = full_B_w(1+diff_L:end-diff_R);
  full_B_e         = full_B_e(1+diff_L:end-diff_R);
  
  
  [full_h_w, ...
   full_h_e ]     = get_edge_values_quadratic ( full_h_P, full_x_P, full_x_w, ...
                                                full_x_e, full_dx_P, ...
                                                full_dx_w, full_dx_e );
                                            
  [full_S_w, ...
   full_S_e ]     = get_edge_values_quadratic ( full_S_P, full_x_P, full_x_w, ...
                                                full_x_e, full_dx_P, ...
                                                full_dx_w, full_dx_e );                                                
                                            
  [ full_dS_dx_w, ...
    full_dS_dx_e ] = get_gradient_values( full_S_P, full_x_P, full_dx_P );    
           


% figure
% plot(full_S_P_R(1:cutoff_R),'r')
% hold on
% plot(full_S_P_L(1:cutoff_L),'g')



% -------------------------------------------------------------------------
% Regrid for a nonuniform grid for transient calculations
% -------------------------------------------------------------------------

 
 % coarse grid:
   coarse_extent_L = 1*full_x_w(1);
   coarse_extent_R = 1*full_x_e(end);
    
   full_x_coarse = [ coarse_extent_L: coarse_dx: coarse_extent_R ] ;
                          
   
 % % fine grid near the termini:   -- not sure this works right.
 % full_x_fine_L = [ full_x_w(1): fine_dx: coarse_extent_L ];
 % full_x_fine_R = [ coarse_extent_R: fine_dx: full_x_e(end) ];
 
 % full_x_edges2 = [ full_x_fine_L full_x_coarse full_x_fine_R ];
   full_x_edges2 = full_x_coarse;
 

  full_x_w2     = full_x_edges2(1:end-1);
  full_x_e2     = full_x_edges2(2:end);
  
  full_dx_w2    = [ (full_x_w2(2) - full_x_w2(1)) diff(full_x_w2)]; 
  full_dx_e2    = [ diff(full_x_e2) (full_x_e2(end) - full_x_e2(end-1)) ];
 
  full_x_P2     = full_x_w2 + full_dx_w2/2;
  
  full_dx_P2    = full_x_e2 - full_x_w2;
  
  
% interpolate everything on to new grid:
% better would be to interpolate everything from x_edges...
  full_W_P         = interp1(full_x_P, full_W_P, full_x_P2, 'linear', 'extrap');
  full_W_w         = interp1(full_x_w, full_W_w, full_x_w2, 'linear', 'extrap');
  full_W_e         = interp1(full_x_e, full_W_e, full_x_e2, 'linear', 'extrap');
  full_slip_P      = interp1(full_x_P, full_slip_P, full_x_P2, 'linear', 'extrap');
  full_slip_w      = interp1(full_x_w, full_slip_w, full_x_w2, 'linear', 'extrap');
  full_slip_e      = interp1(full_x_e, full_slip_e, full_x_e2, 'linear', 'extrap');
  full_A_eff_edges = interp1([full_x_w(1) full_x_e], full_A_eff_edges, [full_x_w2(1) full_x_e2], 'linear', 'extrap');
  full_b_dot_P     = interp1(full_x_P, full_b_dot_P, full_x_P2, 'linear', 'extrap');
  full_b_dot_w     = interp1(full_x_w, full_b_dot_w, full_x_w2, 'linear', 'extrap');
  full_b_dot_e     = interp1(full_x_e, full_b_dot_e, full_x_e2, 'linear', 'extrap');
  full_B_P         = interp1(full_x_P, full_B_P, full_x_P2, 'linear', 'extrap');
  full_B_w         = interp1(full_x_w, full_B_w, full_x_w2, 'linear', 'extrap');
  full_B_e         = interp1(full_x_e, full_B_e, full_x_e2, 'linear', 'extrap');
  full_h_P         = interp1(full_x_P, full_h_P, full_x_P2, 'linear', 'extrap');
  full_h_w         = interp1(full_x_w, full_h_w, full_x_w2, 'linear', 'extrap');
  full_h_e         = interp1(full_x_e, full_h_e, full_x_e2, 'linear', 'extrap');
  full_S_P         = interp1(full_x_P, full_S_P, full_x_P2, 'linear', 'extrap');
  full_S_w         = interp1(full_x_w, full_S_w, full_x_w2, 'linear', 'extrap');
  full_S_e         = interp1(full_x_e, full_S_e, full_x_e2, 'linear', 'extrap');
  full_dS_dx_w     = interp1(full_x_w, full_dS_dx_w, full_x_w2, 'linear', 'extrap');
  full_dS_dx_e     = interp1(full_x_e, full_dS_dx_e, full_x_e2, 'linear', 'extrap');
  
             
  
  % RESET VALUES:
  full_x_P  = full_x_P2;
  full_x_w  = full_x_w2;
  full_x_e  = full_x_e2;
  full_dx_P = full_dx_P2;
  full_dx_w = full_dx_w2;
  full_dx_e = full_dx_e2;
  
  index_div_full = find(full_x_P == 0);   
  full_x_edges   = [ full_x_w(1) full_x_e ]; 

  
  
  

% -------------------------------------------------------------------------
% 4. Replace terminus with a wedge shape
% -------------------------------------------------------------------------


 terminus_cutoff_L = 1;                    % first grid point on left
 terminus_cutoff_R = length(full_x_P);     % last grid point on right

   

 % calculate flux into wedge
 % =========================
   N_use         = 1000;  % number of grid points for flux calculation
   
   x_edges_L     = [full_x_w(1): abs(full_x_w(1))/(N_use-1): 0];
   x_edges_R     = [0: full_x_e(end)/(N_use-1): full_x_e(end)];
   dx_P_L        = diff(x_edges_L);   % new grid
   dx_P_R        = diff(x_edges_R);   % new grid
  
   W_edges_L     = interp1(full_x_edges, [full_W_w(1) full_W_e], x_edges_L);
   W_edges_R     = interp1(full_x_edges, [full_W_w(1) full_W_e], x_edges_R);
   S_dot_L       = zeros(size(x_edges_L(1:end-1)));
   S_dot_R       = zeros(size(x_edges_R(1:end-1)));
   b_dot_edges_L = interp1(full_x_edges, [full_b_dot_w(1) full_b_dot_e], x_edges_L);
   b_dot_edges_R = interp1(full_x_edges, [full_b_dot_w(1) full_b_dot_e], x_edges_R);
   Q_in_L        = 0;
   Q_in_R        = 0; 
   
      
  [ flux_kin_L ] = calc_flux_kin (abs( fliplr(x_edges_L) ), ...
                                  abs( fliplr(x_edges_L) ), ...
                                  fliplr(dx_P_L), fliplr(W_edges_L), ...
                                  fliplr([S_dot_L S_dot_L(end)]), fliplr(b_dot_edges_L), ...
                                  Q_in_L );
                                             
   flux_kin_L     = -fliplr(flux_kin_L);   % flux is negative on the left side                                   
                                  
    
 [ flux_kin_R ] = calc_flux_kin ( x_edges_R, x_edges_R , ...
                                    dx_P_R, W_edges_R, [S_dot_R S_dot_R(end)], ...
                                    b_dot_edges_R, Q_in_R );  



                                
   full_flux_edges_kin = [ flux_kin_L(1:end-1) flux_kin_R ];
   
   full_x_edges_kin    = [ x_edges_L(1:end-1) x_edges_R ];
 
   
 % flux values:  
   Q_in_wedge_L = interp1(full_x_edges_kin, full_flux_edges_kin, ...
                          full_x_e(1));
                      
   Q_in_wedge_R = interp1(full_x_edges_kin, full_flux_edges_kin, ...
                          full_x_w(end));   
 
 
 % use the kinematic value because extrapolated values near the edges are
 % unreliable from the dynamic calculation.
 
 
 
 
% % compare to a dynamic calculation:                                                  
%   full_flux_edges_dyn = calc_flux_dyn( [full_x_w(1) full_x_e], ...
%                                        [full_h_w(1) full_h_e], ...
%                                        [full_dS_dx_w(1) full_dS_dx_e], ...
%                                        [full_slip_w(1) full_slip_e], ...
%                                        [full_W_w(1) full_W_e], ...
%                                        full_A_eff_edges );
 
% % flux values:
%   Q_in_wedge_L = interp1(full_x_edges, full_flux_edges_dyn, ...
%                          full_x_e(1));
%                      
%   Q_in_wedge_R = interp1(full_x_edges, full_flux_edges_dyn, ...
%                          full_x_w(end));
                                        

                              
                                        
  
 % replace ends with a wedge terminus -- calculate ablation rate!
 % ===================================
 [ full_x_P, full_x_w, ...
   full_x_e, full_dx_P, ...
   full_dx_w, full_dx_e, ...
   full_h_P, full_h_w, full_h_e, ...
   full_S_P, full_S_w, full_S_e, ...
   wedge_length_L, ...
   wedge_length_R, ...
   ablation_rate_L, ...
   ablation_rate_R ] = wedge_terminus( full_x_P, full_x_w, full_x_e, ...
                                       full_dx_P, full_dx_w, full_dx_e, ...
                                       full_h_P, full_h_w, full_h_e, ...
                                       full_S_P, full_S_w, full_S_e, ...
                                       full_B_w, full_B_e, ...
                                       Q_in_wedge_L, Q_in_wedge_R, ...
                                       terminus_cutoff_L, ...
                                       terminus_cutoff_R );


% %  prescribe ablation rate!
% %  =========================
% %  
% % Instead of calculate the ablation rate, specify it to be the same
% % as the Paterson model
%   wedge_bdot_L = full_b_dot_w(terminus_cutoff_L); 
%   wedge_bdot_R = full_b_dot_e(terminus_cutoff_R);
% 
%  [ full_x_P, full_x_w, ...
%    full_x_e, full_dx_P, ...
%    full_dx_w, full_dx_e, ...
%    full_h_P, full_h_w, full_h_e, ...
%    full_S_P, full_S_w, full_S_e, ...
%    wedge_length_L, ...
%    wedge_length_R, ...
%    ablation_rate_L, ...
%    ablation_rate_R ] = wedge_terminus2( full_x_P, full_x_w, full_x_e, ...
%                                        full_dx_P, full_dx_w, full_dx_e, ...
%                                        full_h_P, full_h_w, full_h_e, ...
%                                        full_S_P, full_S_w, full_S_e, ...
%                                        full_B_w, full_B_e, ...
%                                        Q_in_wedge_L, Q_in_wedge_R, ...
%                                        terminus_cutoff_L, ...
%                                        terminus_cutoff_R, ...
%                                        wedge_bdot_L, wedge_bdot_R);                                   
  
                                   
                                   
 if ( (wedge_length_L > 10000) || (wedge_length_R > 10000) ) 
     disp( ' ')
     disp(' Wedge terminus length > 10 km -- double check in ldm_extension.m ')
     disp( ' ')
 end
                                   
                                   
% redefine with mass-balance pattern with the wedge:
% ==================================================
  full_b_dot_P = [ ablation_rate_L ...
                   full_b_dot_P(terminus_cutoff_L+1: terminus_cutoff_R-1) ...
                   ablation_rate_R ];
               
  full_b_dot_w = [ ablation_rate_L ...
                   full_b_dot_w(terminus_cutoff_L+1: terminus_cutoff_R) ];
  
  full_b_dot_e = [ full_b_dot_e(terminus_cutoff_L: terminus_cutoff_R-1) ...
                   ablation_rate_R ];
                                     

               



% Check that it is embedded correctly:
plot([full_x_w(1) full_x_e], [full_S_w(1) full_S_e], 'b')
hold on
plot([full_x_w(1) full_x_e], [full_S_w(1) full_S_e], 'b.')
plot([full_x_w(1) full_x_e], [full_B_w(1) full_B_e], 'r')
plot(x_P-x_P(index_ldm_div), h_ldm_P+B_P,'c', 'linewidth', 3)
      
  
               
  

