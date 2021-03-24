% EGR Plotting

cd '/home/users/mborrus/Matlab_HPC'
load('./data/EGR/AM4/AM4_EGR.mat')
mkdir plots/EGR/Comparison

AM4_Data_Path = '/oak/stanford/schools/ees/aditis2/Globus_data/gfdl.intel18-prod-openmp-extra/';
File_Numbers=[2,4];
temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyT_U_SH_O.nc');
lat_AM4 = ncread(temp_path,'grid_yt');
days = (1:100);

load('/scratch/users/mborrus/dycore/h0_EGR_N2/dycore_EGR_N2_00.mat')
egr_mean_h0 = egr_mean;
load('/scratch/users/mborrus/dycore/h4000_EGR_N2/h4000_EGR_N2_00.mat')
egr_mean_h4000 = egr_mean;
load('./data/axis_stuff_64.mat')
lat_dy = lat;
N_Mid_dy = find(lat_dy > 40 & lat_dy < 70);
S_Mid_dy = find(lat_dy < -40 & lat_dy > -70);

N_Mid_am4 = find(lat_AM4 > 40 & lat_AM4 < 70);
S_Mid_am4 = find(lat_AM4 < -40 & lat_AM4 > -70);

%%%%%%%%%

    AM4_dry = squeeze(mean(DRY_egr_mean(1,N_Mid_am4,:),2))*3600*24;
    AM4_moist = squeeze(mean(MOIST_egr_mean(1,N_Mid_am4,:),2))*3600*24;

    h0_dry = squeeze(mean(egr_mean_h0(N_Mid_dy,days),1))*3600*24;
    h4000_dry = squeeze(mean(egr_mean_h4000(N_Mid_dy,days),1))*3600*24;
    
North = figure(1);
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

    AM4_dry = squeeze(DRY_egr_mean(1,:,50))*3600*24;
    AM4_moist = squeeze(MOIST_egr_mean(1,:,50))*3600*24;
    
    h0_dry = squeeze(egr_mean_h0(:,50))*3600*24;
    h4000_dry = squeeze(egr_mean_h4000(:,50))*3600*24;
    
Lat50 = figure(2);
clf, hold on
plot(lat_AM4, AM4_dry)
plot(lat_AM4, AM4_moist)
plot(lat_dy, h0_dry)
plot(lat_dy, h4000_dry)
axis([1 100 0 2.5])
title('50 Lat - EGR Comparison')
legend("AM4 Dry","AM4 Moist", "h0", "h4000")
xlabel("latitude")
ylabel("EGR - 1/days")


saveas(Lat50,['./plots/EGR/Comparison/Lat50_Summary.png'])