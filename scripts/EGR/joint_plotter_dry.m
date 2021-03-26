% EGR Plotting

cd '/home/users/mborrus/Matlab_HPC'
mkdir plots/EGR/Comparison

load('./data/EGR/AM4/AM4_EGR.mat'); clear theta
    lat_am4 = lat; clear lat; p_am4 = p; clear p;
    dtheta_z_dry = squeeze(dtheta_z(2,:,:,:)); clear dtheta_z
    dtheta_p_dry = squeeze(dtheta_p(2,:,:,:)); clear dtheta_p
    N2_dry = squeeze(DRY_N2(2,:,:,:)); clear DRY_N2
    egr_dry = squeeze(DRY_egr(2,:,:,:)); clear DRY_egr
    du_z_am4 = squeeze(du_z(2,:,:,:)); clear du_z
    dtheta_dp_eff = squeeze(dtheta_dp_eff(2,:,:,:)); clear dtheta_dp_eff
    N2_moist = squeeze(MOIST_N2(2,:,:,:)); clear MOIST_N2
    egr_moist = squeeze(MOIST_egr(2,:,:,:)); clear MOIST_egr
    
    N_Mid_am4 = find(lat_am4 > 40 & lat_am4 < 70);
    S_Mid_am4 = find(lat_am4 < -40 & lat_am4 > -70);

load('./data/EGR/dycore/h0_EGR_N2.mat')
    
    egr_h0 = squeeze(egr(2,:,:,:)); 
    N2_h0 = squeeze(N2(2,:,:,:));
    dtheta_z_h0 = squeeze(dtheta_z(2,:,:,:));
    dtheta_p_h0 = squeeze(dtheta_p(2,:,:,:));
    du_z_h0 = squeeze(du_z(2,:,:,:));

load('./data/EGR/dycore/h4000_EGR_N2.mat'); clear theta
    lat_dy = squeeze(lat); clear lat; p_dy = squeeze(p); clear p; 
    egr_h4000 = squeeze(egr(2,:,:,:)); clear egr
    N2_h4000 = squeeze(N2(2,:,:,:)); clear N2
    dtheta_z_h4000 = squeeze(dtheta_z(2,:,:,:)); clear dtheta_z
    dtheta_p_h4000 = squeeze(dtheta_p(2,:,:,:)); clear dtheta_p
    du_z_h4000 = squeeze(du_z(2,:,:,:)); clear du_z

    N_Mid_dy = find(lat_dy > 40 & lat_dy < 70);
    S_Mid_dy = find(lat_dy < -40 & lat_dy > -70);

save('./data/EGR/Comparison_Values.mat')

%%%%%%%%%
%%

% X = lat_dy, Y = p_dy, Input = dtheta_p_h0(:,:,50)
% X_new = lat_am4, Y_new = p_am4

dtheta_p_h0_regrid = interp2(lat_dy,p_dy',dtheta_p_h0(:,:,50)',lat_am4,p_am4');

dtheta_dp_ratio = figure(13);
clf
[M,c] = contour(lat_am4,p_am4, dtheta_p_dry(:,:,50)'./dtheta_p_h0_regrid,'ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dp - am4/h0')
c.LineWidth = 2;
saveas(dtheta_dp_ratio,['./plots/EGR/Comparison/dtheta_dp_ratio.png'])
set(gca, 'YDir','reverse')
!rclone copy /home/users/mborrus/Matlab_HPC/plots remote:plot_folder
%%
%%

N2_plots = figure(12);
clf

subplot(3,1,1)
[M,c] = contour(lat_dy,p_dy, N2_h0(:,:,50)'.*10^5,'ShowText','on');
ylabel('pressure')
title('N2 * 10^5 - h0')
c.LineWidth = 2;
set(gca, 'YDir','reverse')

subplot(3,1,2)
[M,c] = contour(lat_am4',p_am4, N2_dry(:,:,50)'.*10^5,'ShowText','on');
ylabel('pressure')
title('N2 * 10^5 - am4')
c.LineWidth = 2;
set(gca, 'YDir','reverse')

subplot(3,1,3)
N2_h0_regrid = interp2(lat_dy,p_dy',N2_h0(:,:,50)',lat_am4,p_am4');
[M,c] = contour(   lat_am4,p_am4, (N2_dry(:,:,50)'./N2_h0_regrid ),'ShowText','on','LevelList',major);
xlabel('latitude')
ylabel('pressure')
title('N2 - am4/h0')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(N2_plots,['./plots/EGR/Comparison/N2_plots.png'])


%%

dtheta_dz_plot_h4000 = figure(11);
clf
[M,c] = contour(lat_dy,p_dy, dtheta_z_h4000(:,:,50)','ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dz - h4000')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(dtheta_dz_plot_h4000,['./plots/EGR/Comparison/dtheta_dz_h4000.png'])

dtheta_dz_plot_h0 = figure(10);
clf
[M,c] = contour(lat_dy',p_dy, dtheta_z_h0(:,:,50)','ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dz - h0')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(dtheta_dz_plot_h0,['./plots/EGR/Comparison/dtheta_dz_h0.png'])

dtheta_dz_plot_am4 = figure(9);
clf
[M,c] = contour(lat_am4,p_am4, dtheta_z_dry(:,:,50)','ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dz - am4')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(dtheta_dz_plot_am4,['./plots/EGR/Comparison/dtheta_dz_am4.png'])

%%

dtheta_dp_plot_h4000 = figure(8);
clf
[M,c] = contour(lat_dy,p_dy, dtheta_p_h4000(:,:,50)','ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dp - h4000')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(dtheta_dp_plot_h4000,['./plots/EGR/Comparison/dtheta_dp_h4000.png'])

dtheta_dp_plot_h0 = figure(7);
clf
[M,c] = contour(lat_dy',p_dy, dtheta_p_h0(:,:,50)','ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dp - h0')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(dtheta_dp_plot_h0,['./plots/EGR/Comparison/dtheta_dp_h0.png'])

dtheta_dp_plot_am4 = figure(6);
clf
[M,c] = contour(lat_am4,p_am4, dtheta_p_dry(:,:,50)','ShowText','on');
xlabel('latitude')
ylabel('pressure')
title('dtheta/dp - am4')
c.LineWidth = 2;
set(gca, 'YDir','reverse')
saveas(dtheta_dp_plot_am4,['./plots/EGR/Comparison/dtheta_dp_am4.png'])
!rclone copy /home/users/mborrus/Matlab_HPC/plots remote:plot_folder
%%%%%%%%%
%%
CompositePlot_day_50 = figure(5); 
clf
% EGR plot 1/s
subplot(3,1,1)
hold on
title("day 50 - EGR")
ylabel("1/s")
    AM4_dry = squeeze(egr_mean_dry(:,50));
%    AM4_moist = squeeze(egr_mean_moist(:,50));
    h0_dry = squeeze(egr_mean_h0(:,50));
    h4000_dry = squeeze(egr_mean_h4000(:,50));  
plot(lat_am4, AM4_dry)
%plot(lat_am4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 0 3e-05])
legend("AM4 Dry", "h0", "h4000")

subplot(3,1,2)
hold on
title("du/dz")
ylabel("1/s")
    AM4_dry = squeeze(nanmean(du_z_am4(:,:,50),2));
    %AM4_moist = squeeze(nanmean(du_z_am4(:,:,50),2));
    h0_dry = squeeze(nanmean(du_z_h0(:,:,50),2));
    h4000_dry = squeeze(nanmean(du_z_h4000(:,:,50),2));  
plot(lat_am4, AM4_dry)
%plot(lat_am4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 -.001 .005])
%legend("AM4 Dry","AM4 Moist", "h0", "h4000")
legend("AM4 Dry", "h0", "h4000")

subplot(3,1,3)
hold on
title("sqrt(N^2)")
ylabel("1/s")
xlabel("latitude")
    AM4_dry = sqrt(squeeze(N2_mean_dry(:,50)));
%    AM4_moist = sqrt(squeeze(N2_mean_moist(:,50)));
    h0_dry = sqrt(squeeze(N2_mean_h0(:,50)));
    h4000_dry = sqrt(squeeze(N2_mean_h4000(:,50)));  
plot(lat_am4, AM4_dry)
%plot(lat_am4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 0 inf])
legend("AM4 Dry", "h0", "h4000")

saveas(CompositePlot_day_50,['./plots/EGR/Comparison/CompositePlot_day_50.png'])

%%
CompositePlot__moist_day_50 = figure(4); 
clf
% EGR plot 1/s
subplot(3,1,1)
hold on
title("day 50 - EGR")
ylabel("1/s")
    AM4_dry = squeeze(egr_mean_dry(:,50));
    AM4_moist = squeeze(egr_mean_moist(:,50));
    h0_dry = squeeze(egr_mean_h0(:,50));
    h4000_dry = squeeze(egr_mean_h4000(:,50));  
plot(lat_am4, AM4_dry)
plot(lat_am4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 0 3e-05])
legend("AM4 Dry","AM4 Moist", "h0", "h4000")

subplot(3,1,2)
hold on
title("du/dz")
ylabel("1/s")
    AM4_dry = squeeze(nanmean(du_z_am4(:,:,50),2));
    AM4_moist = squeeze(nanmean(du_z_am4(:,:,50),2));
    h0_dry = squeeze(nanmean(du_z_h0(:,:,50),2));
    h4000_dry = squeeze(nanmean(du_z_h4000(:,:,50),2));  
plot(lat_am4, AM4_dry)
plot(lat_am4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 -.001 .005])
legend("AM4 Dry","AM4 Moist", "h0", "h4000")

subplot(3,1,3)
hold on
title("sqrt(N^2)")
ylabel("1/s")
xlabel("latitude")
    AM4_dry = sqrt(squeeze(N2_mean_dry(:,50)));
    %AM4_moist = sqrt(squeeze(
    %N2_mean_moist(:,50)));
    h0_dry = sqrt(squeeze(N2_mean_h0(:,50)));
    h4000_dry = sqrt(squeeze(N2_mean_h4000(:,50)));  
plot(lat_am4, AM4_dry)
plot(lat_am4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 0 inf])
legend("AM4 Dry","AM4 Moist", "h0", "h4000")

saveas(CompositePlot__moist_day_50,['./plots/EGR/Comparison/CompositePlot_moist_day_50.png'])
%%%%%%%%%

    AM4_dry = squeeze(mean(DRY_egr_mean(1,N_Mid_am4,:),2))*3600*24;
    AM4_moist = squeeze(mean(MOIST_egr_mean(1,N_Mid_am4,:),2))*3600*24;

    h0_dry = squeeze(mean(egr_mean_h0(N_Mid_dy,days),1))*3600*24;
    h4000_dry = squeeze(mean(egr_mean_h4000(N_Mid_dy,days),1))*3600*24;
    
North = figure(3);
clf, hold on
plot(days, AM4_dry)
plot(days, AM4_moist)
plot(days, h0_dry)
plot(days, h4000_dry)
axis([1 100 0 2.5])
title('40N-70N - EGR Comparison')
legend("AM4 Dry","AM4 Moist", "h0", "h4000")
xlabel("days")
ylabel("1/days")


saveas(North,['./plots/EGR/Comparison/North_Summary.png'])
%%%%%%%%%

    AM4_dry = squeeze(mean(DRY_egr_mean(1,S_Mid_am4,:),2))*3600*24;
    AM4_moist = squeeze(mean(MOIST_egr_mean(1,S_Mid_am4,:),2))*3600*24;

    h0_dry = squeeze(nanmean(egr_mean_h0(S_Mid_dy,days),1))*3600*24;
    h4000_dry = squeeze(nanmean(egr_mean_h4000(S_Mid_dy,days),1))*3600*24;
    
South = figure(2);
clf, hold on
plot(days, AM4_dry)
plot(days, AM4_moist)
plot(days, h0_dry)
plot(days, h4000_dry)
axis([1 100 0 2.5])
title('40S-70S - EGR Comparison')
legend("AM4 Dry","AM4 Moist", "h0", "h4000")
xlabel("days")
ylabel("1/days")

saveas(South,['./plots/EGR/Comparison/South_Summary.png'])
%%%%%%% day 50

    AM4_dry = squeeze(egr_mean_dry(:,50))*3600*24;
    AM4_moist = squeeze(egr_mean_moist(:,50))*3600*24;
    
    h0_dry = squeeze(egr_mean_h0(:,50))*3600*24;
    h4000_dry = squeeze(egr_mean_h4000(:,50))*3600*24;
    
Lat50 = figure(1);
clf, hold on
plot(lat_AM4, AM4_dry)
plot(lat_AM4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([-90 90 0 2.5])
title('50 Lat - EGR Comparison')
legend("AM4 Dry","AM4 Moist", "h0", "h4000")
xlabel("latitude")
ylabel("EGR - 1/days")


saveas(Lat50,['./plots/EGR/Comparison/Lat50_Summary.png'])