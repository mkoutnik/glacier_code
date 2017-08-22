function [ B_P, B_w, B_e, ...
           dB_dx_w, dB_dx_e, ...
           dB_dx_P, S_modern, ...
           value_to_zero_bed ] = load_bed( x_P, x_w, x_e, dx_P, dx_w, dx_e )


%--------------------------------------------------------------
%
%  set up steady bed elevations B_P(x) at positions x_P
%  where surface S(x,t) is defined
%  interpolate bed to interval midpoints where fluxes are
%  calculated
%
%---------------------------------------------------------------

global DIRECTORY_data

global lower_resolution



addpath(DIRECTORY_data)


% % Option flat bed:
% % ================      
%  elevation = 0;
%  B_P = elevation * ones( size(x_P) );
  

% % Option linear ramp:
% % ===================
% B_P      = interp1( [x_P(1) x_P(end)], [0 1000], x_P );
% S_modern = NaN;


% Values from QGIS along centerline:
% ==================================
 load DH_surf_bed.mat   % Used Mette's data, extracted with profile tool in QGIS
 
B_P_temp          = interp1(Hat_centerline_distance, Hat_bed, x_P);
value_to_zero_bed = abs(min(B_P_temp));

% B_P_zero          = B_P_temp + value_to_zero_bed; %use this if you want
% S_modern = interp1(Hat_centerline_distance, Hat_modern_surface, x_P) + value_to_zero_bed; %use  this if you want the lowest elevation to be zero
% 
% 
% if (lower_resolution == 0)
%   B_P = smooth(B_P_zero, 2)';
%  else
%   B_P = smooth(B_P_zero, 1)';
% end



% Or, use actual values (including negatives, should be ok..?)
B_P               = smooth(B_P_temp,1)';
S_modern = interp1(Hat_centerline_distance, Hat_modern_surface, x_P);



% figure
% plot(x_P/1000, B_P, 'r')
% hold on
% plot(Hat_centerline_distance/1000, Hat_bed,'k')
% plot(x_P/1000, S_modern)
% plot(Hat_centerline_distance/1000, Hat_modern_surface,'k')




% % Values from min search
% % ======================
% disp('Using best min search bed!')
% load best_min_search_bed.mat
% B_P = interp1(x_P_best, B_P_best, x_P, 'linear', 'extrap');


% % Old bit for smoothing bed...
% % ==================================
%  B_P_rough  = interp1(beardmore_x, beardmore_bed_elev, x_P);
%  ws = warning('off','all');  % Turn off warning
%  coeff      = polyfit(x_P, B_P_rough, 19); %16);
%  B_P_smooth = polyval(coeff, x_P);
%  warning(ws)
%  B_P_smooth(end-5:end) = B_P_smooth(end-5);
%  B_P_smooth(1:5)       = B_P_smooth(5); 
%  B_P = B_P_smooth + abs(min(B_P_smooth));   % set lowest value at 0.
 
%   S_rough  = interp1(beardmore_x, beardmore_surface_elev, x_P);
%   ws = warning('off','all');  % Turn off warning
%   coeff    = polyfit(x_P, S_rough, 12); %13);
%   S_smooth = polyval(coeff, x_P);
%   S_modern = S_smooth + abs(min(B_P_smooth));
%   warning(ws)
% figure
% plot(x_P, B_P_rough,'r')
% hold on
% plot(x_P, B_P_smooth,'g')
% plot(x_P, S_rough,'b')
% plot(x_P, S_smooth,'g')

 


% need to extrapolate to control-volume edges and interpolate for
% the interior points
   [B_w, B_e ] = get_edge_values_quadratic ...
                                ( B_P, x_P, x_w, x_e, dx_P, dx_w, dx_e );

%  for ii = 1:length(B_w)
%      if (B_w(ii) < 0 )
%          B_w(ii) = 0;
%      end
%      if (B_e(ii) < 0)
%          B_e(ii) = 0;
%      end
%  end

% the routine get_edge_values_quadratic is at least 10 times faster 
% than interp1 and using 'extrap'


% Bed gradients
% -------------
  [ dB_dx_w, dB_dx_e ] = get_gradient_values( B_P, x_P, dx_P );
       
  
   dB_dx_P = interp1([x_w(1) x_e], [dB_dx_w(1) dB_dx_e], x_P);
   
   