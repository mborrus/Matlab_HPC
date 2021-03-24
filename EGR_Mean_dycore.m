%% Eady growth rate
% Created 11/26/2018 for no EEI paper rewrite

%% Load data
cd '/home/users/mborrus/Matlab_HPC'
load('./data/axis_stuff_64.mat','lat','lon', 'P')
mkdir data/EGR/dycore

OutputNumber = ["0","10"];
DycoreRun = ["h0","h4000"];
DycoreRun_Choice=2;

dycore_Data_Path = '/scratch/users/mborrus/dycore/';

days = (1:100);
Pressure_range = find(P<950 & P>450);
p = P(Pressure_range); 

H = 7300; % scale height, m
z = -H*log(p/1000);
r = 6.37e6; % Radius of Earth
omega = 2*pi/(24*3600);
f = 2*omega*sin(2*pi*lat/360);
g = 9.8;


for run_N = 1:length(OutputNumber)
T = load(strcat(dycore_Data_Path,DycoreRun(DycoreRun_Choice),'/ary',OutputNumber(run_N),'/T_interp_01.mat'),'T_interp_01'); 
u = load(strcat(dycore_Data_Path,DycoreRun(DycoreRun_Choice),'/ary',OutputNumber(run_N),'/u_interp_01.mat'),'u_interp_01'); 
%%
%cut off the top 6 layers of atmosphere
T  = T.T_interp_01(:,:,34:40,:); %pressure from 450hb to 950hb
u  = u.u_interp_01(:,:,34:40,:);
%% Average across pressure and longitude
T = squeeze(mean(T(:,:,:,:),1)); %average across longitude

u = squeeze(mean(u(:,:,:,:),1)); %average across longitude
%%
[la pr d] = size(u);
%% Calculations

% Calculate theta
% we want theta to be the same shape as our other variables (u, v)
% potential temp = T (Reff Pressure/ Pressure)^kappa 
% kappa = 2/7 ratio of gas constant to specific heat capacity
pressure_term = (1000./p).^(2/7);
theta = zeros(la, pr, d);
for i = 1:pr
    theta(:,i,:) = T(:,i,:).*pressure_term(i);
end
%%

% Derivatives WRT z
dtheta_z = zeros(la, pr, d);
du_z = zeros(la, pr, d);
%%
for i=2:pr-1
    dtheta_z(:,i,:) = (theta(:,i+1,:)-theta(:,i-1,:))/...
        (z(i+1)-z(i-1));
    du_z(:,i,:) = (u(:,i+1,:)-u(:,i-1,:))/...
        (z(i+1)-z(i-1));
end
%%
dtheta_z(:,1,:) = (theta(:,2,:)-theta(:,1,:))/(z(2)-z(1));
dtheta_z(:,pr,:) = (theta(:,pr,:)-theta(:,pr-1,:))/(z(pr)-z(pr-1));
du_z(:,1,:) = (u(:,2,:)-u(:,1,:))/(z(2)-z(1));
du_z(:,pr,:) = (u(:,pr,:)-u(:,pr-1,:))/(z(pr)-z(pr-1));
%%
% Buoyancy frequency
% N^2 = g/theta * d/dz(theta)

if exist('N2') == 0
        N2 = zeros(length(OutputNumber),la,pr,d);
        egr = zeros(length(OutputNumber),la,pr,d);
        "vars created"
    else
        "vars already exist"
end

for i=1:pr
    N2(run_N,:,i,:) = (g./theta(:,i,:)) .* dtheta_z(:,i,:);
end
%%
% Eady growth rate
for i = 1:la
    egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(i,:,:))./sqrt(squeeze(N2(run_N,i,:,:))));
end

%%
end

egr_mean = squeeze(nanmean(egr,3));
N2_mean = squeeze(nanmean(N2,3));

save(strcat('./data/EGR/dycore/',DycoreRun(DycoreRun_Choice),'_EGR_N2.mat'), 'egr_mean','N2_mean','dtheta_z','du_z')



