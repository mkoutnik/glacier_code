% Figures to look at results
% ---------------------------

addpath('export_fig')


% Initial:
load save_initial_guess.mat
u_obs           = measures_flowspeed;
u_obs_on_xedges = interp1(measures_centerline_distance, measures_flowspeed, x_edges,'linear', 'extrap'); 
u_calc_ini      = (5/4)*average_vel_estimate;
S_obs           = S_modern;
S_calc_ini      = S_P(1,:);

% Flux:
load save_min_search_flux.mat      % Within expected range
%load save_min_search_flux_fit.mat  % Out of expected range  
u_calc_flux = (5/4)*average_vel_estimate;
S_calc_flux = S_P(1,:);

% Enhancement:
load save_min_search_E.mat         % Within expected range
%load save_min_search_E_fit.mat     % Out of expected range
u_calc_E = (5/4)*average_vel_estimate;
S_calc_E = S_P(1,:);

factor_P_best = factor_P;
x_P_best = x_P;
save best_min_search_E.mat factor_P_best x_P_best

% Bed:
load save_min_search_bed.mat       % Within expected range
%load save_min_search_bed_fit.mat   % Out of expected range
u_calc_bed = (5/4)*average_vel_estimate;
S_calc_bed = S_P(1,:);

B_P_best = B_P;
x_P_best = x_P;
save best_min_search_bed.mat B_P_best x_P_best

% Width:
load save_min_search_width.mat     % Within expected range
% load save_min_search_width_fit.mat % Out of expected range
u_calc_width = (5/4)*average_vel_estimate;
S_calc_width = S_P(1,:);

W_P_best = W_P;
x_P_best = x_P;
save best_min_search_width.mat W_P_best x_P_best

% % sliding factor:
% load save_min_search_fs.mat
% u_calc_fs = (5/4) * average_vel_estimate;
% S_calc_fs = S_P(1,:);
% 


% 
% % Figure 1
% % --------
% % Scatterplot of Surface velocity observed vs. Surface velocity calculated
% figure
% set(gcf, 'Units', 'centimeters','position', [35 20 9 9])
% %subplot('position', [0.25 0.75 0.7 0.2])
% hold on
% plot(u_obs_on_xedges, u_calc_ini,'k.')
% %plot(u_obs_on_xedges, u_calc_flux,'r.')
% plot(u_obs_on_xedges, u_calc_E, 'g.')
% plot(u_obs_on_xedges, u_calc_bed,'b.')
% plot(u_obs_on_xedges, u_calc_width,'m.')
% plot([0 400], [0 400], 'color', [0.6 0.6 0.6])
% %plot(u_obs_on_xedges, u_calc_fs, 'c.')
% xlabel('Observed surface velocity (m/yr)')
% ylabel('Calculated surface velocity (m/yr)')
% 
% % Figure 2
% % --------
% % Scatterplot of Surface observed vs. Surface calculated
% figure(2)
% set(gcf, 'Units', 'centimeters','position', [35 20 9 9])
% %subplot('position', [0.25 0.75 0.7 0.2])
% hold on
% plot(S_obs, S_calc_ini,'k.')
% %plot(S_obs, S_calc_flux, 'r.')
% plot(S_obs, S_calc_E, 'g.')
% plot(S_obs, S_calc_bed, 'b.')
% plot(S_obs, S_calc_width, 'm.')
% %plot(S_obs, S_calc_fs,'c.')
% plot([2300 5500], [2300 5500], 'color', [0.6 0.6 0.6])
% legend('Initial guess', 'Deformation factor', 'Bed', 'Width','location', 'southeast')
% xlabel('Observed relative surface elevation (m)')
% ylabel('Calculated relative surface elevation (m)')
% xlim([2300 5500])
% ylim([2300 5500])



% Figure 1
% --------
% Scatterplot of Surface velocity observed vs. Surface velocity calculated
figure
set(gcf, 'Units', 'centimeters','position', [35 20 12 9])
hold on
plot(x_edges/1000, u_calc_ini-u_obs_on_xedges,'k')
plot(x_edges/1000, u_calc_E-u_obs_on_xedges, 'g')
plot(x_edges/1000, u_calc_bed-u_obs_on_xedges,'c--')
plot(x_edges/1000, u_calc_width-u_obs_on_xedges,'m')
plot([x_edges(1) x_edges(end)]/1000, [0 0], 'color', [0.6 0.6 0.6])
%plot(u_obs_on_xedges, u_calc_fs, 'c.')
xlabel('Distance along flowband (km)', 'fontsize', 14)
ylabel({'(Calculated - Observed)'; 'surface velocity (m/yr)'}, 'fontsize', 14)
legend('Initial guess', '\Delta Deformation', '\Delta Bed', '\Delta Width','location', 'northwest')
ylim([-60 210])
title('Surface velocity mismatch', 'fontweight', 'bold', 'fontsize', 14)

export_fig surface_velocity_mismatch -pdf -transparent


% Figure 2
% --------
% Scatterplot of Surface observed vs. Surface calculated
figure(2)
set(gcf, 'Units', 'centimeters','position', [35 20 12 9])
hold on
plot(x_P/1000, S_calc_ini-S_obs,'k')
plot(x_P/1000, S_calc_E-S_obs, 'g')
plot(x_P/1000, S_calc_bed-S_obs, 'c--')
plot(x_P/1000, S_calc_width-S_obs, 'm')
%plot(S_obs, S_calc_fs,'c.')
plot([x_edges(1) x_edges(end)]/1000, [0 0], 'color', [0.6 0.6 0.6])
legend('Initial guess', '\Delta Deformation', '\Delta Bed', '\Delta Width','location', 'northwest')
xlabel('Distance along flowband (km)', 'fontsize', 14)
ylabel({'(Calculated - Observed)'; 'surface elevation (m)'}, 'fontsize', 14)
title('Surface elevation mismatch', 'fontweight', 'bold', 'fontsize', 14)

export_fig surface_elevation_mismatch -pdf -transparent



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% Plot transient runs
% 1) initial guess
% 2) Best deformation
% 3) Best bed
% 4) Best width


load save_transient_initial_guess.mat
S_P_ini = S_P;

load save_transient_best_bed.mat
S_P_best_bed = S_P;

load save_transient_best_width.mat
S_P_best_width = S_P;

load save_transient_best_E.mat
S_P_best_E = S_P;



figure(3)
set(gcf, 'Units', 'centimeters','position', [5 20 20 10])
% lower left corner coords (x,y), width in pixels, height in pixels

%first_val = 46  % At 14 ka   % =8 is where stops adjusting to initial guess
first_val = N_t_mesh

% a) Grounding line elevation history
subplot('position', [0.1 0.15 0.25 0.75])
gl_pos = 1;
plot(t_P/1000, S_P_ini(:,gl_pos)-S_P_ini(first_val,gl_pos),'k')
%hold on
%plot(t_P/1000, S_P_best_E(:,gl_pos),'g')
%plot(t_P/1000, S_P_best_bed(:,gl_pos),'c--')
%plot(t_P/1000, S_P_best_width(:,gl_pos),'m')
set(gca,'FontSize', 14, 'XTICK', [-20:5:0], 'XTickLabel', {'20','15','10','5','0'})
ylabel({'Surface-elevation change (m)'}, 'fontsize', 14)
xlabel('Time (kyr BP)')
xlim([-14 0])
ylim([0 800])
title('Grounding line', 'fontweight', 'bold', 'fontsize', 14)


% b) Mid-glacier elevation history
subplot('position', [0.4 0.15 0.25 0.75])
mid_pos = round(length(x_P)/2);
plot(t_P/1000, S_P_ini(:,mid_pos)-S_P_ini(first_val,mid_pos),'k')
hold on
plot(t_P/1000, S_P_best_E(:,mid_pos)-S_P_best_E(first_val,mid_pos),'g')
plot(t_P/1000, S_P_best_bed(:,mid_pos)-S_P_best_bed(first_val,mid_pos),'c--')
plot(t_P/1000, S_P_best_width(:,mid_pos)-S_P_best_width(first_val,mid_pos),'m')
set(gca,'FontSize', 14, 'YTICK', [0:200:800], 'YTickLabel', {'','','','','','','','',''},'XTICK', [-20:5:0], 'XTickLabel', {'20','15','10','5','0'})
xlabel('Time (kyr BP)')
title('Middle of glacier', 'fontweight', 'bold', 'fontsize', 14)
xlim([-14 0])
ylim([0 800])

% c) Head of glacier elevation history
subplot('position', [0.7 0.15 0.25 0.75])
head_pos = length(x_P)-2;
plot(t_P/1000, S_P_ini(:,head_pos)-S_P_ini(first_val,head_pos),'k')
hold on
plot(t_P/1000, S_P_best_E(:,head_pos)-S_P_best_E(first_val,head_pos),'g')
plot(t_P/1000, S_P_best_bed(:,head_pos)-S_P_best_bed(first_val,head_pos),'c--')
plot(t_P/1000, S_P_best_width(:,head_pos)-S_P_best_width(first_val,head_pos),'m')
h = legend('Initial guess', '\Delta Deformation', '\Delta Bed', '\Delta Width','location', 'northwest')
set(h, 'fontsize', 12)
set(gca,'FontSize', 14, 'YTICK', [-800:200:0], 'YTickLabel', {'','','','','','','','',''},'XTICK', [-20:5:0], 'XTickLabel', {'20','15','10','5','0'})
xlabel('Time (kyr BP)')
title('Head of glacier', 'fontweight', 'bold', 'fontsize', 14)
xlim([-14 0])
ylim([0 800])


export_fig surface_elevation_histories -pdf -transparent





