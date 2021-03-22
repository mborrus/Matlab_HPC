%% Eady growth rate
% Created 11/26/2018 for no EEI paper rewrite

%% Load data
clear
load('./Data/axis_stuff_64.mat','lat','lon', 'P')
OutputNumber = '00';
T = load(['/Volumes/SD/h4000/',OutputNumber,'/T_interp_01.mat'],'T_interp_01'); 
u = load(['/Volumes/SD/h4000/',OutputNumber,'/u_interp_01.mat'],'u_interp_01'); 
%%
%cut off the top 6 layers of atmosphere
T  = T.T_interp_01(:,:,34:40,:); %pressure from 450hb to 950hb
u  = u.u_interp_01(:,:,34:40,:);
%% Average across pressure and longitude
T = squeeze(mean(T(:,:,:,:),1)); %average across longitude

u = squeeze(mean(u(:,:,:,:),1)); %average across longitude
%%
p = P(34:40); 
clear P

H = 7300; % scale height, m
z = -H*log(p/1000);
r = 6.37e6; % Radius of Earth
omega = 2*pi/(24*3600);
f = 2*omega*sin(2*pi*lat/360);
g = 9.8;

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
N2 = zeros(la,pr,d);
for i=1:pr
    N2(:,i,:) = (g./theta(:,i,:)) .* dtheta_z(:,i,:);
end
%%
% Eady growth rate
egr = zeros(la,pr,d);
for i = 1:la
    egr(i,:,:) = abs(f(i)*du_z(i,:,:)./sqrt(N2(i,:,:)));
end
%%
egr_mean = squeeze(mean(egr,2));

N2_mean = squeeze(mean(N2,2));

%%
save(['h4000_EGR_N2_',OutputNumber,'.mat'], 'egr_mean','N2_mean')
%%
%EGR plot across time
% cmap = jet(15);
% names = [];
% figure(1),clf,hold on
% for i = 0:8
%     plot(lat,egr_mean(:,1+49*i),'Color', cmap(i+1, :))
%     names = [names ; 1+49*i];
% end
% legend(num2str(names(:)))
% title('EGR over time - 20')
% xlabel('lat')
% ylabel('egr')
% set(gca, 'fontsize', 14)
% 
% 
% %%
% %EGR plots at +-60deg (-11 and +54)
% figure(2),clf,hold on
% subplot(3,1,1)
% plot((1:400),egr_mean_lon_pr(11,:),'b')
% title('egr at -60 deg lat - 20')
% ylabel('egr')
% set(gca, 'fontsize', 14)
% 
% subplot(3,1,2)
% plot((1:400),egr_mean_lon_pr(54,:), 'r')
% title('egr at +60 deg lat')
% ylabel('egr')
% set(gca, 'fontsize', 14)
% 
% subplot(3,1,3), hold on
% title('comparison')
% plot((1:400),egr_mean_lon_pr(54,:), 'r')
% plot((1:400),egr_mean_lon_pr(11,:),'b')
% legend('+60','-60')
% xlabel('days')
% ylabel('egr')
% set(gca, 'fontsize', 14)
% 





