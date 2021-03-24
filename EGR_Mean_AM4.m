%% Eady growth rate for AM4 - Moist and Dry
% Script to be run in sherlock to make plots
% Marshall Borrus
cd '/home/users/mborrus/Matlab_HPC'
AM4_Data_Path = '/oak/stanford/schools/ees/aditis2/Globus_data/gfdl.intel18-prod-openmp-extra/';
mkdir data/EGR/AM4

%%
prompt = 'Do you want plots? 1 = yes, 0 = no';
x = input(prompt)
ploton = x; 
% Get the pressure and lats values
File_Numbers=[2,4];

temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyT_U_SH_O.nc');
pressure_levels = ncread(temp_path,'pfull');
lat = ncread(temp_path,'grid_yt');
days = (1:100);
Pressure_range = find(pressure_levels<950 & pressure_levels>450);
p = pressure_levels(Pressure_range); 

H = 7300; % scale height, m
z = -H*log(p/1000);
r = 6.37e6; % Radius of Earth
omega = 2*pi/(24*3600);
f = 2*omega*sin(2*pi*lat/360);
g = 9.8;
lambda =.6;

for run_N = 1:length(File_Numbers)

    temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(run_N)),'/dailyT_U_SH_O.nc');
    
    %Load U and T
    start_loc = [ 1 1 Pressure_range(1) 1];
    count = [inf inf length(Pressure_range) 100];
    u = ncread(temp_path,'ucomp',start_loc,count);
    T = ncread(temp_path,'temp',start_loc,count);
    
    %Zonal mean
    T = squeeze(mean(T(:,:,:,:),1)); %average across longitude
    u = squeeze(mean(u(:,:,:,:),1)); %average across longitude
    [la pr d] = size(u);
    
    pressure_term = (1000./p).^(2/7);
    theta = zeros(la, pr, d);
    for j = 1:pr
        theta(:,j,:) = T(:,j,:).*pressure_term(j);
    end
    
    dtheta_z = zeros(la, pr, d);
    du_z = zeros(la, pr, d);

    for i=2:pr-1
    dtheta_z(:,i,:) = (theta(:,i+1,:)-theta(:,i-1,:))/...
        (z(i+1)-z(i-1));
    du_z(:,i,:) = (u(:,i+1,:)-u(:,i-1,:))/...
        (z(i+1)-z(i-1));
    end
    
    dtheta_z(:,1,:) = (theta(:,2,:)-theta(:,1,:))/(z(2)-z(1));
    dtheta_z(:,pr,:) = (theta(:,pr,:)-theta(:,pr-1,:))/(z(pr)-z(pr-1));
    du_z(:,1,:) = (u(:,2,:)-u(:,1,:))/(z(2)-z(1));
    du_z(:,pr,:) = (u(:,pr,:)-u(:,pr-1,:))/(z(pr)-z(pr-1));
    
    if exist('DRY_N2') == 0;
        DRY_N2 = zeros(length(File_Numbers),la,pr,d);
        DRY_egr = zeros(length(File_Numbers),la,pr,d);
        MOIST_N2 = zeros(length(File_Numbers),la,pr,d);
        MOIST_egr = zeros(length(File_Numbers),la,pr,d);
        "vars created"
    else
        "vars already exist"
    end
    
    for i=1:pr
        DRY_N2(run_N,:,i,:) = (g./theta(:,i,:)) .* dtheta_z(:,i,:);
    end
    
    for i = 1:la
        DRY_egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(i,:,:))./sqrt(squeeze(DRY_N2(run_N,i,:,:))));
    end
    
    dtheta_dp_eff = zeros(la,pr,d);
    for lat_N = 1:la
        for day_N = 1:d
            dtheta_dp_eff(lat_N,:,day_N) = eff_stat_stab(p', T(lat_N,:,day_N), lambda);
        end
    end
    
    for i=1:pr
        MOIST_N2(run_N,:,i,:) = (g./theta(:,i,:)) .* dtheta_dp_eff(:,i,:);
    end
    
    for i = 1:la
        MOIST_egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(i,:,:))./sqrt(squeeze(MOIST_N2(run_N,i,:,:))));
    end
    
end

    DRY_egr_mean = squeeze(mean(DRY_egr,3));
    DRY_N2_mean = squeeze(mean(DRY_N2,3));

    MOIST_egr_mean = squeeze(mean(MOIST_egr,3));
    MOIST_N2_mean = squeeze(mean(MOIST_N2,3));
    
    DRY_dtheta_z = dtheta_z;
    MOIST_dtheta_z = dtheta_dp_eff;
    save(['./data/EGR/AM4/AM4_EGR.mat'], 'DRY_egr_mean','DRY_N2_mean','DRY_dtheta_z','du_z','MOIST_dtheta_z','MOIST_egr_mean','MOIST_N2_mean')
    
    
    
    %% 
    if ploton
    ratio_dry_over_moist = DRY_egr_mean./MOIST_egr_mean;

       
    Global_average = mean(ratio_dry_over_moist(2,:,:),2);
    global_avg = figure(1);
    plot(1:100,squeeze(Global_average))
    title('Global Average Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    axis([ 1 100 0 7])
    saveas(global_avg,['./plots/EGR/AM4/Global_AVG.png'])
    
    sixty = figure(2);
    
    plot(1:100,squeeze(ratio_dry_over_moist(2,75,:)))
    title('60N Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    axis([ 1 100 0 7])
    saveas(sixty,['./plots/EGR/AM4/60N.png'])
    
    Comparison = figure(3);
    hold on;
    plot(1:100,squeeze(ratio_dry_over_moist(2,75,:)))
    plot(1:100,squeeze(Global_average))
    title('Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    legend('60N','Global')
    axis([ 1 100 0 7])
    saveas(Comparison,['./plots/EGR/AM4/Comparison.png'])
    
    
%%
    Global_average = mean(ratio_dry_over_moist(:,:,:),2); hold on;
    global_avg = figure(1);
    plot(1:100,squeeze(Global_average))
    title('Global Average Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    legend('Global (1)','Global (2)')
    axis([ 1 100 0 10])
    saveas(global_avg,['./plots/EGR/AM4/Global_AVG_2.png'])
    
    Mid_Lats = figure(2); hold on;
    lat_ranges = lat;
    N_Mid = find(lat_ranges > 20 & lat_ranges < 50);
    S_Mid = find(lat_ranges < -20 & lat_ranges > -50);
    plot(1:100,squeeze(mean(ratio_dry_over_moist(:,N_Mid,:),2)))
    plot(1:100,squeeze(mean(ratio_dry_over_moist(:,S_Mid,:),2)))
    title('Mid-Lats Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    axis([ 1 100 0 10])
    legend('20-50N (1)','20-50N (2)','20-50S (1)','20-50S (2)')
    saveas(Mid_Lats,['./plots/EGR/AM4/Mid_Lats_20_50.png'])
    
    Mid_Lats_2 = figure(20); hold on;
    lat_ranges = lat;
    N_Mid = find(lat_ranges > 40 & lat_ranges < 70);
    S_Mid = find(lat_ranges < -40 & lat_ranges > -70);
    plot(1:100,squeeze(mean(ratio_dry_over_moist(:,N_Mid,:),2)))
    plot(1:100,squeeze(mean(ratio_dry_over_moist(:,S_Mid,:),2)))
    title('Mid-Lats Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    axis([ 1 100 0 10])
    legend('40-70N (1)','40-70N (2)','40-70S (1)','40-70S (2)')
    saveas(Mid_Lats_2,['./plots/EGR/AM4/Mid_Lats_40_70.png'])
    
    Comparison = figure(3); hold on;
    hold on;
    plot(1:100,squeeze(mean(ratio_dry_over_moist(:,N_Mid,:),2)))
    plot(1:100,squeeze(Global_average))
    title('Global vs MidLats Dry EGR/Moist EGR')
    xlabel('days')
    ylabel('dry/moist')
    legend('20-50N (1)','20-50N (2)','Global (1)','Global (2)')
    axis([ 1 100 0 10])
    saveas(Comparison,['./plots/EGR/AM4/Comparison_2.png'])
    end

