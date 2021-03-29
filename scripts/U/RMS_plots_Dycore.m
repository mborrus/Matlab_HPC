%% RMS_plots V2 
% Script to be run in sherlock to make plots
% Marshall Borrus
% This is where the data is held

Dycore_Data_Path = '/scratch/users/mborrus/h0/';
cd /home/users/mborrus/Matlab_HPC
mkdir plots/RMS     %Where these specific plots should go

%%

% Get the pressure and lats values
% size(P) = 40; size(lat) = 64; size(lon) = 128;

load(strcat(Dycore_Data_Path,'axis_stuff_64.mat'),'lat','lon','P');
days = (1:100);
    
% Set pressure and lat ranges
whole_P_range = find(P<850 & P>1);
pressure_levels = P(whole_P_range);  clear P;
Nlon = length(lon);
Nlat = length(lat);
Np   = length(pressure_levels);
Ntime= length(days);
Nruns = 21; Nrun = Nruns-1;  %%%%%%%%% CHANGE THIS IF YOU NEED TO %%%%%%%%%%%
lat_ranges = lat;     clear lat;


Lower_Tropo = find(pressure_levels<850 & pressure_levels>500);
Upper_Tropo = find(pressure_levels<500 & pressure_levels>100);
Stratosphere= find(pressure_levels<100 & pressure_levels>1);

N_Hi = find(lat_ranges >  50);
S_Hi = find(lat_ranges < -50);
N_Mid = find(lat_ranges > 20 & lat_ranges < 50);
S_Mid = find(lat_ranges < -20 & lat_ranges > -50);
Equator = find(lat_ranges > -20 & lat_ranges < 20);

%% Grab the U values from each folder

        % ucomp = Size:       128x64x40x400
        % Dimensions: lon, lat, P, days

        % Get initial values so we don't need to save more than we need
    temp_path = strcat(Dycore_Data_Path,'ary',num2str(0),'/u_interp_01.mat');
    temp_u = load(temp_path,'u_interp_01');
    U_0 = temp_u.u_interp_01(:,:,whole_P_range,1:100);
    
    clear temp_path temp_u
    U_diff = zeros(20,Nlon,Nlat,Np,Ntime);
for i = 1:Nrun
            
            % The files are spread out in 20 different folders, this gets you to each of them
    temp_path = strcat(Dycore_Data_Path,'ary',num2str(i),'/u_interp_01.mat');
    temp_u = load(temp_path,'u_interp_01');
    temp_u = temp_u.u_interp_01(:,:,whole_P_range,1:100);
    
            % subtract temp_u
    U_diff(i,:,:,:,:) = temp_u - U_0;
end
    clear U_1 temp_u

N_low_T = length(Lower_Tropo);
N_Up_T = length(Upper_Tropo);
N_Strat = length(Stratosphere);

%%
c = 'Color';
c1 = ["#a1dab4","#41b6c4","#225ea8"];

% GLOBAL
        
        range_temp = 1:Nlat; rng = length(range_temp);
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

    axis([1 100 0 75]), title('Global RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/Global_h0.png'])

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

    axis([1 100 0 75]), title('North Hi-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/North_Hi_h0.png'])
    
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

    axis([1 100 0 75]), title('North Mid-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/North_Mid_h0.png'])

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

    axis([1 100 0 75]), title('Low-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/Low_h0.png'])

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

    axis([1 100 0 75]), title('South Mid-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/South_Mid_h0.png'])

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

    axis([1 100 0 75]), title('South Hi-lat RMS - Dycore h0'), xlabel('Days'), ylabel('RMSE')

    saveas(MainPlot,['./plots/RMS/South_Hi_h0.png'])


%%

