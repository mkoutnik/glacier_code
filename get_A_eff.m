function [ A_eff_P, ...
           A_eff_w, A_eff_e ] = get_A_eff ( T, A_0, ...
                                            zhat_grid, d_zhat, ...
                                            x_P, x_w, x_e, dx_P, ...
                                            dx_w, dx_e )                                   
                                       
                                       
% -------------------------------------------------------------------------
% calculate the effective isothermal softness parameter 

% this code is taken from velocity_field.m

% -------------------------------------------------------------------------


global Q R n
                                            
     

 zhat_grid_use = flipud(zhat_grid);   % now goes surface to bed in z.
 T_use         = flipud(T); 

% integrate over z
% -----------------
  int_exp_QRT = cumtrapz( (exp(-Q ./ (R*T_use))) .* (1-zhat_grid_use).^n ) * d_zhat;  


% integrate over z again
% -----------------------
  int_int_exp_QRT = trapz( int_exp_QRT ) * d_zhat;  % trapz assumes dz=1.

                                       
                                       
                                       
% ----------------------------------------------------
% effective isothermal softness parameter, A_eff(x,t)
% ----------------------------------------------------
  A_eff_P = A_0 .* (n+2) .* int_int_exp_QRT;
                                       
                                       
                                       
% values at edges:     
% ================     
   [ A_eff_w, ...
     A_eff_e ] = get_edge_values_quadratic ( A_eff_P, x_P, x_w, x_e, ...
                                             dx_P, dx_w, dx_e );                                     
                                       
                                       
                                       
                                       
                                       



                                            
                                            
                                            
% old bits, might be junk. (judged this 5 Sept 2009)
                                       
                                       
% function [ A_eff_xt, A_T_z ] = get_A_eff ( x_P, S_P, B_P, W_P, ...
%                                            T_field, K_vec, c_vec, x_therm, ...
%                                            z_therm )
% 
% % ------------------------------------------------------------------------
% % Michelle Koutnik
% % with notes/bits from Tom Neumann and Ed Waddington
% % April 2007
% 
% % parts of code taken from tm_geometry.m
% 
% % from a starting surface temperature profile and ice thickness, 
% % get the temperature at depth using a heat-balance model (Firestone et al.,
% % 1990)
% % Get the effective isothermal softness parameter that gives the same 
% % ice velocity and flux as temperature-dependent softness parameter
% 
% % -------------------------------------------------------------------------
% 
% 
%  global s_per_year  rho_ice  g  n  A_0  Q  R  E  Qgeo
%  global N_t_mesh  N_x_mesh
% 
% 
% 
% A_eff_xt = NaN * ones(N_t_mesh, N_x_mesh);
% 
% 
% for itime = 1:N_t_mesh
% 
% ice_thickness = S_P(itime,:) - B_P;
% bed_elev      = B_P; % output from ss model
% surface_elev  = S_P(itime,:);
% width         = W_P;  % output of ss model
%  
% T_z_K         = T_field(:, :, itime);
% x             = x_therm(:, :, itime);
% z             = z_therm(:, :, itime);
% 
% K             = K_vec(:,itime)/ s_per_year;  % W/mK 
% c             = c_vec(:,itime);
% 
% kappa         = K./(rho_ice*c);	% thermal diffusivity (m^2/s)
% 
%  
% 
% 
% % first get A(T(x,z)), using temperature T_z(z,x)
% % -----------------------------------------------
% 
% A_T_z = A_0*exp(-Q./(R.*T_z_K));
% 
% 
% % % z_hat = (z - B) / (S - B)
% icethick_temp = interp1( x_P, ice_thickness, x(1,:), 'linear', 'extrap');
% ice_thick_grid = repmat( icethick_temp, length(x(:,1)), 1); 
% bed_temp = interp1( x_P, bed_elev, x(1,:), 'linear', 'extrap');
% bed_grid = repmat( bed_temp, length(x(:,1)), 1);
% zhat_temp = (z - bed_grid) ./ (ice_thick_grid);
% 
% zhat_pos = repmat(abs(min(zhat_temp')), length(zhat_temp(1,:)), 1)' + zhat_temp;
% 
% zhat = zhat_pos;
% 
% 
% % % setup non-dimensional z, zhat
% % for j=1:length(y(1,:))
% %   max_thickness(j) = max(y(:,j));
% %   zhat(:,j) = y(:,j) ./ max_thickness(j);
% % 
% %   nan_vals = isnan(A_T_z(:,j));
% %   if (length(nan_vals) > 0)
% %       nan_index = find(nan_vals == 1);
% %       zhat(nan_index,j) = 0;
% %       y(nan_index,j) = 0;
% %       A_T_z(nan_index,j) = 0;
% %   end
%   
%   
% end
% 
% % Tom notes:
%  % z(1) = 0 at bed, positive upward -- seems like z=0 at surface and
%  % positive downward to bed...
% 
%  
% % get effective isothermal softness parameter
% % --------------------------------------------
% 
% % take first integral from 0 to zhat
% % cumtrapz assumes unit spacing, really use dzhat...
% 
% for ii = 1:length(z(1,:))
%     
%  dzhat_temp(:,ii) = diff(zhat(:,ii));
%  
%  index1 = find(dzhat_temp(:,ii) == -1);
%  dzhat_temp(index1,ii) = dzhat_temp(index1-1,ii);
%  
% end
% 
% dzhat = [zeros(1,length(zhat(1,:)))' dzhat_temp']';
% 
% 
% % multiply by dzhat because cumtrapz assumes unit spacing.
% F = cumtrapz( A_T_z .* ( (1-zhat).^n )) .* dzhat;
% 
% 
% % now get ready for second integral over zhat, 0-1
% F_end = NaN*ones(1,length(z(1,:)));
% for k = 1:length(z(1,:))
%     index_end = find(z(:,k) == max(z(:,k)));
%     F_end(k) = F(index_end, k);
% end
% 
% % F(b) - F(a)
% int_A_T_z = (F_end - F(1,:));
% 
% 
% 
% % at last!
% % ---------
% 
% A_eff = (n+2) * int_A_T_z;
% 
% 
% % convert to Pa^-3 yr^-1
% A_eff = A_eff * (s_per_year / 1e9 );
% 
% 
% 
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % compare to linear profile:
% % 
% % T_surf = T_s(1,1);
% % 
% %       dz  = 1;                         %  m
% %       z_T  = 800;                     %  m
% %       z = ( 0: dz: z_T );
% %       n_z = length(z);
% %       T = NaN*ones(length(z), 1);
% %       T(1) = T_surf;
% % 
% %       for i_z = 2:n_z
% %           T(i_z) = T(i_z-1) + ( (Qgeo/K)  )*dz;
% %           
% %       end    %  for i_z = 1:n_z ...
% %       
% %   
% %   %  find softness A(T) at bottom of profile
% %       
% %       A_T = A_0 * exp( -Q/(R * (T(end)+273.15) ) );  % kPa^-n s^-1
% %    
% %    % convert to Pa^-3 mars_year^-1
% %       A_T = A_T * (sec_per_year/1e09);
% %       
% 
% 
% 
% 
