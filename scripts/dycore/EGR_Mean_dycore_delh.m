%% Eady growth rate
% Created 11/26/2018 for no EEI paper rewrite

%% Load data
delh = ["10", "20", "30", "40", "50", "60", "70"];
pert = ["0", "1", "2"];

cd '/home/users/mborrus/Matlab_HPC'
load('./data/axis_stuff_64.mat','lat','lon', 'P')
mkdir data/EGR/dycore

lat = lat'; 
%pert = ["0","10"];
%DycoreRun = ["h0","h4000"];
%DycoreRun_Choice=2;

dycore_Data_Path = '/scratch/users/mborrus/dycore/';

days = (1:100);
Pressure_range = find(P<950 & P>200);
p = P(Pressure_range); 

H = 7300; % scale height, m
z = -H*log(p/1000);
r = 6.37e6; % Radius of Earth
omega = 2*pi/(24*3600);
f = 2*omega*sin(2*pi*lat/360);
g = 9.8;

%for delh_N = 3:length(delh)
for delh_N = 2
    clear theta dtheta_z dtheta_p du_z N2 egr
for pert_N = 1:length(pert)
    
T = load(strcat(dycore_Data_Path,delh(delh_N),'/',pert(pert_N),'/T_interp_01.mat'),'T_interp_01'); 
u = load(strcat(dycore_Data_Path,delh(delh_N),'/',pert(pert_N),'/u_interp_01.mat'),'u_interp_01'); 
%%
%cut off the top 6 layers of atmosphere
T  = T.T_interp_01(:,:,Pressure_range,:); %pressure from 450hb to 950hb
u  = u.u_interp_01(:,:,Pressure_range,:);
%% Average across pressure and longitude
T = squeeze(mean(T(:,:,:,:),1)); %average across longitude
u = squeeze(mean(u(:,:,:,:),1)); %average across longitude
%%
[la pr d] = size(u);
if exist('theta') == 0;
        theta = zeros(length(pert),la, pr, d);
        dtheta_z = zeros(length(pert),la, pr, d);
        dtheta_p = zeros(length(pert),la, pr, d);
        du_z = zeros(length(pert),la, pr, d);
        N2 = zeros(length(pert),la,pr,d);
        egr = zeros(length(pert),la,pr,d);
            "vars created"
    else
            "vars already exist"
end

%% Calculations

% Calculate theta
% we want theta to be the same shape as our other variables (u, v)
% potential temp = T (Reff Pressure/ Pressure)^kappa 
% kappa = 2/7 ratio of gas constant to specific heat capacity
pressure_term = (1000./p).^(2/7);
    for j = 1:pr
        theta(pert_N,:,j,:) = T(:,j,:).*pressure_term(j);       
    end
            "theta ran"


    for i=2:pr-1
    dtheta_z(pert_N,:,i,:) = (theta(pert_N,:,i+1,:)-theta(pert_N,:,i-1,:))/...
        (z(i+1)-z(i-1));
    du_z(pert_N,:,i,:) = (u(:,i+1,:)-u(:,i-1,:))/...
        (z(i+1)-z(i-1));
    dtheta_p(pert_N,:,i,:) = (theta(pert_N,:,i+1,:)-theta(pert_N,:,i-1,:))/...
        (p(i+1)-p(i-1));         
    end
            "derivatives ran"
    
    dtheta_p(pert_N,:,1,:) = (theta(pert_N,:,2,:)-theta(pert_N,:,1,:))/(p(2)-p(1));
    dtheta_p(pert_N,:,pr,:) = (theta(pert_N,:,pr,:)-theta(pert_N,:,pr-1,:))/(p(pr)-p(pr-1));
    dtheta_z(pert_N,:,1,:) = (theta(pert_N,:,2,:)-theta(pert_N,:,1,:))/(z(2)-z(1));
    dtheta_z(pert_N,:,pr,:) = (theta(pert_N,:,pr,:)-theta(pert_N,:,pr-1,:))/(z(pr)-z(pr-1));
    du_z(pert_N,:,1,:) = (u(:,2,:)-u(:,1,:))/(z(2)-z(1));
    du_z(pert_N,:,pr,:) = (u(:,pr,:)-u(:,pr-1,:))/(z(pr)-z(pr-1));
%%
% Buoyancy frequency
% N^2 = g/theta * d/dz(theta)



for i=1:pr
    N2(pert_N,:,i,:) = (g./theta(pert_N,:,i,:)) .* dtheta_z(pert_N,:,i,:);
end
%%
% Eady growth rate
for i = 1:la
    egr(pert_N,i,:,:) = abs(f(i)*squeeze(du_z(pert_N,i,:,:))./sqrt(squeeze(N2(pert_N,i,:,:))));
end
delh(delh_N)
pert(pert_N)
%%
end

% egr_mean = squeeze(nanmean(egr,3));
% N2_mean = squeeze(nanmean(N2,3));

% save(strcat('./data/EGR/dycore/',DycoreRun(DycoreRun_Choice),'_EGR_N2.mat'), 'egr_mean','N2_mean','dtheta_z','du_z')
save(strcat('./data/EGR/dycore/',delh(delh_N),'_EGR_N2.mat'), 'lat', 'p', 'theta','dtheta_z','dtheta_p','du_z','N2', 'egr');
    "run saved"
end


