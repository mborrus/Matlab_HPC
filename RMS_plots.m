%% RMS_plots V2 
% Script to be run in sherlock to make plots
% Marshall Borrus
% This is where the data is held

AM4_Data_Path = '/scratch/users/mborrus/Globus_data/gfdl.intel18-prod-openmp-extra/';
mkdir plots         %where plots should go
mkdir plots/RMS     %Where these specific plots should go

%%

% Get the pressure and lats values

temp_path = strcat(AM4_Data_Path,num2str(1),'/dailyU.nc');
pressure_levels = ncread(temp_path,'pfull');
lat_ranges = ncread(temp_path,'grid_yt');
days = (1:100);

% Set pressure and lat ranges

Lower_Tropo = find(pressure_levels<850 & pressure_levels>500);
Upper_Tropo = find(pressure_levels<500 & pressure_levels>100);
Stratosphere= find(pressure_levels<100 & pressure_levels>1);



N_Hi = find(lat_ranges >  50);
S_Hi = find(lat_ranges < -50);
N_Mid = find(lat_ranges > 20 & lat_ranges < 50);
S_Mid = find(lat_ranges < -20 & lat_ranges > -50);
Equator = find(lat_ranges > -20 & lat_ranges < 20);

%% Grab the U values from each folder

        % ucomp = Size:       144x90x33x365
        % Dimensions: grid_xt,grid_yt,pfull,time
        % Dimensions: lon, lat,pfull,time

        % Get initial values so we don't need to save more than we need
    temp_path = strcat(AM4_Data_Path,num2str(1),'/dailyU.nc');
    temp_u = ncread(temp_path,'ucomp');
    
    U_1 = temp_u(:,:,1:24,1:100);
    clear temp_path temp_u
    U_diff = zeros(19,144,90,24,100);
for i = 2:20
            
            % The files are spread out in 20 different folders, this gets you to each of them
    temp_path = strcat(AM4_Data_Path,num2str(i),'/dailyU.nc');
    temp_u = ncread(temp_path,'ucomp');
    
            % subtract temp_u
    U_diff(i-1,:,:,:,:) = temp_u(:,:,1:24,1:100) - U_1;
end
    clear U_1 temp_u

[Nrun,Nlon,Nlat,Np,Ntime] = size(U_diff);
N_low_T = length(Lower_Tropo);
N_Up_T = length(Upper_Tropo);
N_Strat = length(Stratosphere);

%%
c = 'Color';
c1 = ["#a1dab4","#41b6c4","#225ea8"];

% GLOBAL
        
        range_temp = 1:90; rng = length(range_temp);
        U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
        U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
        U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
        
    MainPlot = figure(1)
    clf, hold on
    plot_lt = plot(days, U_LT, c, c1(1)); 
    plot_ut = plot(days, U_UT, c, c1(2));
    plot_s = plot(days, U_S , c, c1(3));
    plots = [plot_lt(1), plot_ut(1), plot_s(1)];
    hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

    axis([1 100 0 75]), title('Global RMS'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/Global.png'])

% N_Hi
        clear U_LT U_UT U_S
        range_temp = N_Hi; rng = length(range_temp);
        U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
        U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
        U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
    MainPlot = figure(1)
    clf, hold on
    plot_lt = plot(days, U_LT, c, c1(1)); 
    plot_ut = plot(days, U_UT, c, c1(2));
    plot_s = plot(days, U_S , c, c1(3));
    plots = [plot_lt(1), plot_ut(1), plot_s(1)];
    hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

    axis([1 100 0 75]), title('North Hi-lat RMS'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/North_Hi.png'])
    
% N_Mid
        clear U_LT U_UT U_S
        range_temp = N_Mid; rng = length(range_temp);
        U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
        U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
        U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
    MainPlot = figure(1)
    clf, hold on
    plot_lt = plot(days, U_LT, c, c1(1)); 
    plot_ut = plot(days, U_UT, c, c1(2));
    plot_s = plot(days, U_S , c, c1(3));
    plots = [plot_lt(1), plot_ut(1), plot_s(1)];
    hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

    axis([1 100 0 75]), title('North Mid-lat RMS'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/North_Mid.png'])

% Low_lat
        clear U_LT U_UT U_S
        range_temp = Equator; rng = length(range_temp);
        U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
        U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
        U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
    MainPlot = figure(1)
    clf, hold on
    plot_lt = plot(days, U_LT, c, c1(1)); 
    plot_ut = plot(days, U_UT, c, c1(2));
    plot_s = plot(days, U_S , c, c1(3));
    plots = [plot_lt(1), plot_ut(1), plot_s(1)];
    hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

    axis([1 100 0 75]), title('Low-lat RMS'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/Low.png'])

% Mid_S_lat
        clear U_LT U_UT U_S
        range_temp = S_Mid; rng = length(range_temp);
        U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
        U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
        U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
    MainPlot = figure(1)
    clf, hold on
    plot_lt = plot(days, U_LT, c, c1(1)); 
    plot_ut = plot(days, U_UT, c, c1(2));
    plot_s = plot(days, U_S , c, c1(3));
    plots = [plot_lt(1), plot_ut(1), plot_s(1)];
    hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

    axis([1 100 0 75]), title('South Mid-lat RMS'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/South_Mid.png'])

% High_S_lat
        clear U_LT U_UT U_S
        range_temp = S_Hi; rng = length(range_temp);
        U_LT = squeeze(rms(reshape(U_diff(:,:,range_temp,Lower_Tropo,:), Nrun,rng*Nlon*N_low_T,Ntime),2));
        U_UT = squeeze(rms(reshape(U_diff(:,:,range_temp,Upper_Tropo,:), Nrun,rng*Nlon*N_Up_T,Ntime),2));
        U_S =  squeeze(rms(reshape(U_diff(:,:,range_temp,Stratosphere,:),Nrun,rng*Nlon*N_Strat,Ntime),2));
    MainPlot = figure(1)
    clf, hold on
    plot_lt = plot(days, U_LT, c, c1(1)); 
    plot_ut = plot(days, U_UT, c, c1(2));
    plot_s = plot(days, U_S , c, c1(3));
    plots = [plot_lt(1), plot_ut(1), plot_s(1)];
    hleg = legend(plots,'Lower Troposphere','Upper Troposphere','Stratosphere');

    axis([1 100 0 75]), title('South Hi-lat RMS'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/South_Hi.png'])


%%

