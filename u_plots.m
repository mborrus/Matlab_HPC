%% Script to be run in sherlock to make plots
AM4_Data_Path = '/scratch/users/mborrus/Globus_data/gfdl.intel18-prod-openmp-extra/';
mkdir plots
%% Grab the U values from each folder
% ucomp = Size:       144x90x33x365
%         Dimensions: grid_xt,grid_yt,pfull,time
%         Dimensions: lon, lat,pfull,time
U_850HPa_60lat = []; %pressure level = 23: 848.8 HPa
U_500HPa_60lat = []; %pressure level = 18: 532.5 HPa
U_10HPa_60lat = []; %pressure level = 3: 10.7 HPa

% Lat(75) = 59 deg
% 
for i = 1:20;
    temp_path = strcat(AM4_Data_Path,num2str(i),'/dailyU.nc');
    temp_u = ncread(temp_path,'ucomp');
    %
    U_850HPa_60lat(i,:,:) = squeeze(temp_u(:,75,23,:));
    U_500HPa_60lat(i,:,:) = squeeze(temp_u(:,75,18,:));
    U_10HPa_60lat(i,:,:) = squeeze(temp_u(:,75,3,:));
end

days = (1:365);
%%
Single850 = figure(1);
clf, hold on
plot(days,squeeze(U_850HPa_60lat(:,1,:)))
axis([1 100 -inf inf])
title('U at 1.25 lon - 60 deg N - 850 HPa')
xlabel('Days')
ylabel('u, m/s')

saveas(Single850,['./plots/single_850.png'])
%%
All850 = figure(2);
clf, hold on
data = mean(U_850HPa_60lat(:,:,:),2);
plot(days,squeeze(data(:,:)))
axis([1 100 -inf inf])
title('Zonal Mean')
xlabel('Days')
ylabel('Zonal Mean u - 60 deg N - 850 HPa')

saveas(All850,['./plots/all_850.png'])
%%
Single500 = figure(10);
clf, hold on
plot(days,squeeze(U_500HPa_60lat(:,1,:)))
axis([1 100 -inf inf])
title('U at 1.25 lon - 60 deg N - 500 HPa')
xlabel('Days')
ylabel('u, m/s')


saveas(Single500,['./plots/single_500.png'])
%%
All500 = figure(20);
clf, hold on
data = mean(U_500HPa_60lat(:,:,:),2);
plot(days,squeeze(data(:,:)))
axis([1 100 -inf inf])
title('Zonal Mean')
xlabel('Days')
ylabel('Zonal Mean u - 60 deg N - 500 HPa')

saveas(All500,['./plots/all_500.png'])
%%
Single10 = figure(100);
clf, hold on
plot(days,squeeze(U_10HPa_60lat(:,1,:)))
axis([1 100 -inf inf])
title('U at 1.25 lon - 60 deg N - 10 HPa')
xlabel('Days')
ylabel('u, m/s')

saveas(Single10,['./plots/single_10.png'])
%%
All10 = figure(200);
clf, hold on
data = mean(U_10HPa_60lat(:,:,:),2);
plot(days,squeeze(data(:,:)))
axis([1 100 -inf inf])
title('Zonal Mean')
xlabel('Days')
ylabel('Zonal Mean u - 60 deg N - 10 HPa')

saveas(All10,['./plots/all_10.png'])