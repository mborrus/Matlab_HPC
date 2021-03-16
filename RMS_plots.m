%% Script to be run in sherlock to make plots
%This is where the data is held
AM4_Data_Path = '/scratch/users/mborrus/Globus_data/gfdl.intel18-prod-openmp-extra/';
mkdir plots %where plots should go
mkdir plots/RMS %Where these specific plots should go
%% Grab the U values from each folder
%       ucomp = Size:       144x90x33x365
%         Dimensions: grid_xt,grid_yt,pfull,time
%         Dimensions: lon, lat,pfull,time
U_850HPa_60lat = []; %pressure level = 23: 848.8 HPa
U_500HPa_60lat = []; %pressure level = 18: 532.5 HPa
U_10HPa_60lat = []; %pressure level = 3: 10.7 HPa
days = (1:365);
% Lat(75) = 59 deg
% 
for i = 1:20;
    % The files are spread out in 20 different folders, this gets you to each of them
    temp_path = strcat(AM4_Data_Path,num2str(i),'/dailyU.nc');
    temp_u = ncread(temp_path,'ucomp');
    % This creates three seperate i x 144 x 365 variables
    % This could be done in a single variable, but then it would show up in
    % your workspace as a 4-D double (gross!)
    U_850HPa_60lat(i,:,:) = squeeze(temp_u(:,75,23,:));
    U_500HPa_60lat(i,:,:) = squeeze(temp_u(:,75,18,:));
    U_10HPa_60lat(i,:,:) = squeeze(temp_u(:,75,3,:));
end
%% Index Matrix
% Create an index of what to subtract from your base in each iteration, so
% as to not include itself (not sure how much this effects it to start
% with...)
Index(1,:) = [2:20];
Index(20,:) = [1:19];
for i = 2:19
    Index(i,:)=[1:i-1 i+1:20];
end
%%
for i = 1:20
    % This takes error between the base and the 19 alternatives for the
    % first longitude bin, with is 1.25 deg
    error = squeeze((U_850HPa_60lat(i,1,:)-U_850HPa_60lat(Index,1,:)));
    % Gets the RMS of the error
    rmserror = rms(error);
    % Save the RMS error to the corresponding run
    RMS850single(i,:)=rmserror;
    
    %Same as above, but this takes the zonal mean of the u_comp prior to
    %calculating error
    errormean = squeeze((mean(U_850HPa_60lat(i,:,:),2)-mean(U_850HPa_60lat(Index,1,:),2)));
    rmserrormean = rms(errormean);
    RMS850mean(i,:)=rmserrormean;
    
    %repeate for 500 Hpa, etc... 
    error = squeeze((U_500HPa_60lat(i,1,:)-U_500HPa_60lat(Index,1,:)));
    rmserror = rms(error);
    RMS500single(i,:)=rmserror;
    
    errormean = squeeze((mean(U_500HPa_60lat(i,:,:),2)-mean(U_500HPa_60lat(Index,1,:),2)));
    rmserrormean = rms(errormean);
    RMS500mean(i,:)=rmserrormean;
    
    error = squeeze((U_10HPa_60lat(i,1,:)-U_10HPa_60lat(Index,1,:)));
    rmserror = rms(error);
    RMS10single(i,:)=rmserror;
    
    errormean = squeeze((mean(U_10HPa_60lat(i,:,:),2)-mean(U_10HPa_60lat(Index,1,:),2)));
    rmserrormean = rms(errormean);
    RMS10mean(i,:)=rmserrormean;
end
%%
Single850 = figure(1);
clf, hold on
plot(days,RMS850single)
axis([1 100 -inf inf])
title('RMSE at 1.25 lon - 60 deg N - 850 HPa')
xlabel('Days')
ylabel('RMSE')

saveas(Single850,['./plots/RMS/single_850.png'])
%%
All850 = figure(2);
clf, hold on
plot(days,RMS850mean)
axis([1 100 -inf inf])
title('Zonal Mean RMSE - 60 deg N - 850 HPa')
xlabel('Days')
ylabel('RMSE')

saveas(All850,['./plots/RMS/all_850.png'])
%%
Single500 = figure(10);
clf, hold on
plot(days,RMS500single)
axis([1 100 -inf inf])
title('RMSE at 1.25 lon - 60 deg N - 500 HPa')
xlabel('Days')
ylabel('RMSE')


saveas(Single500,['./plots/RMS/single_500.png'])
%%
All500 = figure(20);
clf, hold on
plot(days,RMS500mean)
axis([1 100 -inf inf])
title('Zonal Mean RMSE - 60 deg N - 500 HPa')
xlabel('Days')
ylabel('RMSE')

saveas(All500,['./plots/RMS/all_500.png'])
%%
Single10 = figure(100);
clf, hold on
plot(days,RMS10single)
axis([1 100 -inf inf])
title('RMSE at 1.25 lon - 60 deg N - 10 HPa')
xlabel('Days')
ylabel('RMSE')

saveas(Single10,['./plots/RMS/single_10.png'])
%%
All10 = figure(200);
clf, hold on
plot(days,RMS10mean)
axis([1 100 -inf inf])
title('Zonal Mean RMSE - 60 deg N - 10 HPa')
xlabel('Days')
ylabel('RMSE')

saveas(All10,['./plots/RMS/all_10.png'])