
clear
cd '/home/users/mborrus/Matlab_HPC'
addpath('/home/users/mborrus/Matlab_HPC/scripts/EGR')
addpath('/home/users/mborrus/Matlab_HPC/plots')
addpath('/home/users/mborrus/Matlab_HPC/scripts/cbrewer')
load('./data/EGR/lambda.mat')
%%
load('AM4_Data.mat')

% Create the data paths for AM4 data
    AM4_Data.path{1} = '/scratch/users/mborrus/AM4/Base/';
    AM4_Data.path{2} = '/scratch/users/mborrus/AM4/co2/';
    AM4_Data.path{3} = '/scratch/users/mborrus/AM4/quarter/';
    AM4_Data.path{4} = '/scratch/users/mborrus/AM4/m3/';
    AM4_Data.path{5} = '/scratch/users/mborrus/AM4/m2/';
    AM4_Data.path{6} = '/scratch/users/mborrus/AM4/m1/';
    AM4_Data.path{7} = '/scratch/users/mborrus/AM4/p1/';
    AM4_Data.path{8} = '/scratch/users/mborrus/AM4/p2/';
    AM4_Data.path{9} = '/scratch/users/mborrus/AM4/p3/';
    AM4_Data.path{10} = '/scratch/users/mborrus/AM4/SST_p4/';
    AM4_Data.path{11} = '/scratch/users/mborrus/AM4/SST_m4/';
    AM4_Data.path{12} = '/scratch/users/mborrus/AM4/SST_p8/';
    AM4_Data.path{13} = '/scratch/users/mborrus/AM4/SST_p12/';
    AM4_Data.path{14} = '/scratch/users/mborrus/AM4/SST_p16/';


% Create run names, # of files, title names, and plot colors    
    AM4_Data.run_type{1} = '+0C Base';      AM4_Data.File_Numbers{1} = [1:20]; AM4_Data.codename{1} = 'base'; AM4_Data.Color{1}= [0.6980    0.8745    0.5412; 0.2000    0.6275    0.1725]; AM4_Data.Line_Style{1}='-';
    AM4_Data.run_type{2} = '+4C Very Hot';  AM4_Data.File_Numbers{2} = [1:20]; AM4_Data.codename{2} = 'p4'; AM4_Data.Color{2}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{2}='-.';
    AM4_Data.run_type{3} = '-4C Very Cold'; AM4_Data.File_Numbers{3} = [1:20]; AM4_Data.codename{3} = 'm4'; AM4_Data.Color{3}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]; AM4_Data.Line_Style{3}='-.';
    AM4_Data.run_type{4} = '-3C Cold'; AM4_Data.File_Numbers{4} = [1:2]; AM4_Data.codename{4} = 'm3'; AM4_Data.Color{4}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.9; AM4_Data.Line_Style{4}=':';
    AM4_Data.run_type{5} = '-2C Colder'; AM4_Data.File_Numbers{5} = [1:2]; AM4_Data.codename{5} = 'm2'; AM4_Data.Color{5}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.8; AM4_Data.Line_Style{5}=':'; 
    AM4_Data.run_type{6} = '-1C Coldish'; AM4_Data.File_Numbers{6} = [1:2]; AM4_Data.codename{6} = 'm1'; AM4_Data.Color{6}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.7; AM4_Data.Line_Style{6}=':';
    AM4_Data.run_type{7} = '+1C Warmish'; AM4_Data.File_Numbers{7} = [1:2]; AM4_Data.codename{7} = 'p1'; AM4_Data.Color{7}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.9; AM4_Data.Line_Style{7}=':';
    AM4_Data.run_type{8} = '+2C Warmer'; AM4_Data.File_Numbers{8} = [1:2]; AM4_Data.codename{8} = 'p2'; AM4_Data.Color{8}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.8; AM4_Data.Line_Style{8}=':';
    AM4_Data.run_type{9} = '+3C Warm'; AM4_Data.File_Numbers{9} = [1:2]; AM4_Data.codename{9} = 'p3'; AM4_Data.Color{9}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]*.7; AM4_Data.Line_Style{9}=':';
    AM4_Data.run_type{10} = '+4C QOB (hot)'; AM4_Data.File_Numbers{10} = [1:4]; AM4_Data.codename{10} = 'p4 QOB'; AM4_Data.Color{10}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{10}='-';
    AM4_Data.run_type{11} = '-4C QOB (cold)'; AM4_Data.File_Numbers{11} = [1:4]; AM4_Data.codename{11} = 'm4 QOB'; AM4_Data.Color{11}= [0.6510    0.8078    0.8902; 0.1216    0.4706    0.7059]; AM4_Data.Line_Style{11}='-';
    AM4_Data.run_type{12} = '+8C QOB'; AM4_Data.File_Numbers{12} = [1:3]; AM4_Data.codename{12} = 'p8 QOB'; AM4_Data.Color{12}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{12}='-';
    AM4_Data.run_type{13} = '+12C QOB'; AM4_Data.File_Numbers{13} = [1:3]; AM4_Data.codename{13} = 'p12 QOB'; AM4_Data.Color{13}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{13}='-';
    AM4_Data.run_type{14} = '+16C QOB'; AM4_Data.File_Numbers{14} = [1:3]; AM4_Data.codename{14} = 'p16 QOB'; AM4_Data.Color{14}= [0.9843    0.6039    0.6000; 0.8902    0.1020    0.1098]; AM4_Data.Line_Style{14}='-';

% Anxillary data: Lat ranges, pressure ranges 
    %Lat Ranges
    AM4_Data.lat_range_names{1} = "North High"; AM4_Data.lat_range{1} = find(AM4_Data.lat >  50);
    AM4_Data.lat_range_names{2} = "South High"; AM4_Data.lat_range{2} = find(AM4_Data.lat <  -50);
    AM4_Data.lat_range_names{3} = "North Mid"; AM4_Data.lat_range{3} = find(AM4_Data.lat > 20 & AM4_Data.lat < 50);
    AM4_Data.lat_range_names{4} = "South Mid"; AM4_Data.lat_range{4} = find(AM4_Data.lat < -20 & AM4_Data.lat > -50);
    AM4_Data.lat_range_names{5} = "Equator"; AM4_Data.lat_range{5} = find(AM4_Data.lat > -20 & AM4_Data.lat < 20);
    AM4_Data.lat_range_names{6} = "45 deg N"; AM4_Data.lat_range{6} = find(AM4_Data.lat > 42.5 & AM4_Data.lat < 47.5);
    %Pressure Ranges
    AM4_Data.p_range_names{1} = "Lower Tropo"; AM4_Data.p_range{1} = find(AM4_Data.p<850 & AM4_Data.p>500); AM4_Data.p_color{1} = "#a1dab4";
    AM4_Data.p_range_names{2} = "Upper Tropo"; AM4_Data.p_range{2} = find(AM4_Data.p<500 & AM4_Data.p>100); AM4_Data.p_color{2} = "#41b6c4";
    AM4_Data.p_range_names{3} = "Stratosphere"; AM4_Data.p_range{3} = find(AM4_Data.p<100 & AM4_Data.p>1); AM4_Data.p_color{3} = "#225ea8";
    
    save(['AM4_Data.mat'],'AM4_Data')
    
%% 2. Variable Checking
%%%% A. Mean, STD, and RMS Variable Creation
    %U - wind speeds
    for run_selection = [12:14];
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels>1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,24,100);
                var_long = zeros(length(File_Numbers),144,90,24,100);
            temp_var = ncread(temp_path,'ucomp');
            var_1 = temp_var(:,:,1:24,1:100);
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
            temp_var = ncread(temp_path,'ucomp');
            var_diff(i-1,:,:,:,:) = temp_var(:,:,1:24,1:100) - var_1;
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100); i
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff); clear var_RMS var_Mean var_STD
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp); pres
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,p_range{pres},:), Nrun,rng*Nlon*length(p_range{pres}),Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,p_range{pres},:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,p_range{pres},:),1,[1,2,3,4],'omitnan'));
            end
        end
            AM4_Data.U_RMS{run_selection} = var_RMS;
            AM4_Data.U_Mean{run_selection} = var_Mean;
            AM4_Data.U_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% T - Temperatures
    for run_selection = [12:14];
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels>1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,24,100);
                var_long = zeros(length(File_Numbers),144,90,24,100);
            temp_var = ncread(temp_path,'temp');
            var_1 = temp_var(:,:,1:24,1:100);
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
            temp_var = ncread(temp_path,'temp');
            var_diff(i-1,:,:,:,:) = temp_var(:,:,1:24,1:100) - var_1;
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100); i
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff);
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp); pres
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,p_range{pres},:), Nrun,rng*Nlon*length(p_range{pres}),Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,p_range{pres},:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,p_range{pres},:),1,[1,2,3,4],'omitnan'));
            end
        end
            AM4_Data.T_RMS{run_selection} = var_RMS;
            AM4_Data.T_Mean{run_selection} = var_Mean;
            AM4_Data.T_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% P - Precipitation
    for run_selection = [12:14];
        File_Numbers=AM4_Data.File_Numbers{run_selection};
        AM4_Data_Path = AM4_Data.path{run_selection};
        temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        
        % Get the pressure and lats values    
        pressure_levels = ncread(temp_path,'pfull');
            p_range{1} = find(pressure_levels<850 & pressure_levels>500);
            p_range{2} = find(pressure_levels<500 & pressure_levels>100);
            p_range{3} = find(pressure_levels<100 & pressure_levels>1);    
        lat_ranges = ncread(temp_path,'grid_yt');
        for i = 1:length(File_Numbers)
            if i == 1
                var_diff = zeros(length(File_Numbers)-1,144,90,100);
                var_long = zeros(length(File_Numbers),144,90,100);
            temp_var = ncread(temp_path,'precip');
            var_1 = temp_var(:,:,1:100);
            var_long(i,:,:,:) = temp_var(:,:,1:100);
            else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(i)),'/dailyUTP.nc');
            temp_var = ncread(temp_path,'precip');
            var_diff(i-1,:,:,:) = temp_var(:,:,1:100) - var_1;
            var_long(i,:,:,:) = temp_var(:,:,1:100); i
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Ntime] = size(var_diff);
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp); pres
                var_RMS(pres,lats,:) = squeeze(nanmean(rms(reshape(var_diff(:,:,range_temp,:), Nrun,rng*Nlon,Ntime),2),1));
                var_Mean(pres,lats,:) = squeeze(nanmean(var_long(:,:,range_temp,:),[1,2,3,4]));
                var_STD(pres,lats,:) = squeeze(std(var_long(:,:,range_temp,:),1,[1,2,3,4],'omitnan'));
            end
        end
            AM4_Data.P_RMS{run_selection} = var_RMS;
            AM4_Data.P_Mean{run_selection} = var_Mean;
            AM4_Data.P_STD{run_selection} = var_STD;
            run_selection
            clear var_diff var_long var_RMS var_Mean var_STD
    end
    save(['AM4_Data.mat'],'AM4_Data');
    %% Single Value Means
    for run_selection = [12:14];
        for lats = 1:6
            for pres=1:3
                AM4_Data.U_value{run_selection}(pres,lats)=squeeze(nanmean(AM4_Data.U_Mean{run_selection}(pres,lats,:),[1,2,3]));
                AM4_Data.T_value{run_selection}(pres,lats)=squeeze(nanmean(AM4_Data.T_Mean{run_selection}(pres,lats,:),[1,2,3]));
                AM4_Data.P_value{run_selection}(pres,lats)=squeeze(nanmean(AM4_Data.P_Mean{run_selection}(pres,lats,:),[1,2,3]));
            end
        end
    end
    save(['AM4_Data.mat'],'AM4_Data');
%%%% B. Mean+STD and RMS ploting
%% Mean + STD plotting
    %% U - Wind Ppeed Sat Time
    for run_selection = [12:14]; 
        sat_times = zeros(3,6);
        for lat_band = [1:6];
            File_Numbers=AM4_Data.File_Numbers{run_selection};
            savepath = strcat('./plots/Master/RMS/',AM4_Data.codename{run_selection},'/'); mkdir(savepath)
            run_name = AM4_Data.run_type{run_selection};
            range_temp = 1:90; rng = length(range_temp);
            U_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.U_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                plot_rms(i) = plot(1:100, a, 'Color', AM4_Data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
                plot_rms(i) = plot(1:100, b, 'Color', AM4_Data.p_color{i},'LineWidth',2);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
                plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
                i;
            end
            hleg = legend(plot_rms(:),AM4_Data.p_range_names{:}); 
            axis([1 100 0 75]), title([AM4_Data.lat_range_names{lat_band},' RMS :',run_name]), xlabel('Days'), ylabel('RMSE')
            saveas(U_RMS_Plot,strcat(savepath,'_U_',AM4_Data.lat_range_names{lat_band},'.png'))     
            sat_times(:,lat_band)=b_sat(:) ;
            run_selection 
        end
        AM4_Data.U_Sat_Time{run_selection}=sat_times;
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %% T - Temp Sat Time
    for run_selection = [12:14]; 
        sat_times = zeros(3,6);
        for lat_band = [1:6];
            File_Numbers=AM4_Data.File_Numbers{run_selection};
            savepath = strcat('./plots/Master/RMS/',AM4_Data.codename{run_selection},'/'); mkdir(savepath);
            run_name = AM4_Data.run_type{run_selection};
            range_temp = 1:90; rng = length(range_temp);
            U_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.T_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                plot_rms(i) = plot(1:100, a, 'Color', AM4_Data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
                plot_rms(i) = plot(1:100, b, 'Color', AM4_Data.p_color{i},'LineWidth',2);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
                plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
            end
            hleg = legend(plot_rms(:),AM4_Data.p_range_names{:}); 
            axis([1 100 0 75]), title([AM4_Data.lat_range_names{lat_band},' T RMS :',run_name]), xlabel('Days'), ylabel('RMSE')
            saveas(U_RMS_Plot,strcat(savepath,'_T_',AM4_Data.lat_range_names{lat_band},'.png'))     
            sat_times(:,lat_band)=b_sat(:) ;
        end
        run_selection
        AM4_Data.T_Sat_Time{run_selection}=sat_times;
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %% P - Precip Sat Time
    for run_selection = [12:14]; 
        sat_times = zeros(3,6);
        for lat_band = [1:6];
            File_Numbers=AM4_Data.File_Numbers{run_selection};
            savepath = strcat('./plots/Master/RMS/',AM4_Data.codename{run_selection},'/'); mkdir(savepath);
            run_name = AM4_Data.run_type{run_selection};
            range_temp = 1:90; rng = length(range_temp);
            U_RMS_Plot = figure(1); clf, hold on, clear a b plot_rms b_sat  
            for i = 1:3; 
                a =squeeze(AM4_Data.P_RMS{run_selection}(i,lat_band,:)); b = smooth(a,10);
                plot_rms(i) = plot(1:100, a, 'Color', AM4_Data.p_color{i},'LineWidth',2); plot_rms(i).Color(4)=.3;
                plot_rms(i) = plot(1:100, b, 'Color', AM4_Data.p_color{i},'LineWidth',2);
                b_sat(i) = find(diff(b)./b(2:end)<.03,1)+1;
                plot(b_sat(i),b(b_sat(i)),'r*','LineWidth',2)
            end
            hleg = legend(plot_rms(:),AM4_Data.p_range_names{:}); 
            axis([1 100 0 inf]), title([AM4_Data.lat_range_names{lat_band},' P RMS :',run_name]), xlabel('Days'), ylabel('RMSE')
            saveas(U_RMS_Plot,strcat(savepath,'_P_',AM4_Data.lat_range_names{lat_band},'.png'))     
            sat_times(:,lat_band)=b_sat(:) ;
        end
        run_selection
        AM4_Data.P_Sat_Time{run_selection}=sat_times;
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %% Scatter Plots - Delete or Move if not interesting 
        for lat_band = 6
            U_RMS_Sat_Time_Scatter = figure(1); clf; hold on; plot_count = 1; legend_names = [];
            savepath = strcat('./plots/Master/Sat_Scatter/',AM4_Data.codename{run_selection},'/'); mkdir(savepath)
            clear x_data y_data
            for plot_choices = [1:11];
                x_data(plot_count,:)=AM4_Data.T_Mean{plot_count}(:,lat_band)-AM4_Data.T_Mean{1}(:,lat_band);y_data(plot_count,:)=AM4_Data.U_Sat_Time{plot_choices}(:,lat_band);
                plot_count=plot_count+1;
            end
            scatter(x_data(:,1),y_data(:,1),'r')
            %scatter(x_data(:,2),y_data(:,1),'b')
            %scatter(x_data(:,3),y_data(:,1),'g')
            ylabel("days"); xlabel("Temp from base"); set(gca,'FontSize',12)
            title(["Saturation Times", AM4_Data.lat_range_names{lat_band}])
            saveas(U_RMS_Sat_Time_Scatter,strcat(savepath,'_U_',AM4_Data.lat_range_names{lat_band},'.png'))     
        end
%% 3. EGR Calculation
    
    for run_selection = [12:14];
        if run_selection>0
    clear lambda_term T u lat p theta dtheta_dz dtheta_dp du_z DRY_N2 DRY_egr MOIST_N2 MOIST_egr dtheta_dz_eff dtheta_dp_eff
    File_Numbers=AM4_Data.File_Numbers{run_selection};
    AM4_Data_Path = AM4_Data.path{run_selection};

        if run_selection >= 10;
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTP.nc');
        else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(1)),'/dailyUTco2.nc');
        end

    pressure_levels = ncread(temp_path,'pfull');
    lat = ncread(temp_path,'grid_yt');
    days = (1:100);
    Pressure_range = find(pressure_levels<850 & pressure_levels>200);
    p = pressure_levels(Pressure_range); 

    H = 7300; % scale height, m
    z = -H*log(p/1000);
    omega = 2*pi/(24*3600);
    f = 2*omega*sin(2*pi*lat/360);
    g = 9.8;

    lambda_input = interp1(lambda_base(:,1),lambda_base(:,2),lat);

    for run_N = 1:length(File_Numbers)
    %for run_N = 1

        if run_selection >= 10;
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(run_N)),'/dailyUTP.nc');
        else
            temp_path = strcat(AM4_Data_Path,num2str(File_Numbers(run_N)),'/dailyUTco2.nc');
        end

        %Load U and T
        start_loc = [ 1 1 Pressure_range(1) 1];
        count = [inf inf length(Pressure_range) 100];
        u = ncread(temp_path,'ucomp',start_loc,count);
        T = ncread(temp_path,'temp',start_loc,count);

        %Zonal mean
        T = squeeze(mean(T(:,:,:,:),1)); %average across longitude
        u = squeeze(mean(u(:,:,:,:),1)); %average across longitude
        [la pr d] = size(u);

        if exist('theta') == 0;
            theta = zeros(length(File_Numbers),la, pr, d);
            dtheta_dz = zeros(length(File_Numbers),la, pr, d);
            dtheta_dp = zeros(length(File_Numbers),la, pr, d);
            du_z = zeros(length(File_Numbers),la, pr, d);
            DRY_N2 = zeros(length(File_Numbers),la,pr,d);
            DRY_egr = zeros(length(File_Numbers),la,pr,d);
            MOIST_N2 = zeros(length(File_Numbers),la,pr,d);
            MOIST_egr = zeros(length(File_Numbers),la,pr,d);
            dtheta_dz_eff = zeros(length(File_Numbers),la,pr,d);
            dtheta_dp_eff = zeros(length(File_Numbers),la,pr,d);
            lambda_term = zeros(length(File_Numbers),la,pr,d);

                "vars created"
        else
                "vars already exist"
        end

        %DRY RUNS 
        lambda = 0;
        for lat_N = 1:la
            for day_N = 1:d
                [dtheta_dz(run_N,lat_N,:,day_N),dtheta_dp(run_N,lat_N,:,day_N),lambda_term(run_N,lat_N,:,day_N)] = eff_stat_stab(p', T(lat_N,:,day_N), lambda);
            end
        end

        %MOIST RUNS

        for lat_N = 1:la
        lambda = lambda_input(lat_N);
            for day_N = 1:d
                [dtheta_dz_eff(run_N,lat_N,:,day_N),dtheta_dp_eff(run_N,lat_N,:,day_N),lambda_term(run_N,lat_N,:,day_N)] = eff_stat_stab(p', T(lat_N,:,day_N), lambda);
            end
        end

        pressure_term = (1000./p).^(2/7);

        for j = 1:pr
            theta(run_N,:,j,:) = T(:,j,:).*pressure_term(j);       
        end
                "theta ran"

        for i=2:pr-1
        du_z(run_N,:,i,:) = (u(:,i+1,:)-u(:,i-1,:))/...
            (z(i+1)-z(i-1));       
        end
                "derivatives ran"

        du_z(run_N,:,1,:) = (u(:,2,:)-u(:,1,:))/(z(2)-z(1));
        du_z(run_N,:,pr,:) = (u(:,pr,:)-u(:,pr-1,:))/(z(pr)-z(pr-1));


        for i=1:pr
            DRY_N2(run_N,:,i,:) = (g./theta(run_N,:,i,:)) .* dtheta_dz(run_N,:,i,:);
        end

        for i = 1:la
            DRY_egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(run_N,i,:,:))./sqrt(squeeze(DRY_N2(run_N,i,:,:))));
        end

        for i=1:pr
            MOIST_N2(run_N,:,i,:) = (g./theta(run_N,:,i,:)) .* dtheta_dz_eff(run_N,:,i,:);
        end

        for i = 1:la
            MOIST_egr(run_N,i,:,:) = abs(f(i)*squeeze(du_z(run_N,i,:,:))./sqrt(squeeze(MOIST_N2(run_N,i,:,:))));
        end
        run_N
    end

    save([strcat('./data/EGR/AM4/AM4_',AM4_Data.codename{run_selection},'.mat')],'lambda_term', 'T','u', 'lat', 'p', 'theta','dtheta_dz','dtheta_dp','du_z','DRY_N2', 'DRY_egr','MOIST_N2','MOIST_egr','dtheta_dz_eff','dtheta_dp_eff');
    AM4_Data.dry_egr{run_selection} = DRY_egr;
    AM4_Data.moist_egr{run_selection} = MOIST_egr;
    AM4_Data.dry_N2{run_selection} = DRY_N2;
    AM4_Data.moist_N2{run_selection} = MOIST_N2;
    "saved"
        end
    end
    save(['AM4_Data.mat'],'AM4_Data')
    %
    % Mean value
    for run_selection = [12:14];
        if run_selection > 0;
            for i = 1:6
                dry_means(i)=nanmean(AM4_Data.dry_egr{run_selection}(:,AM4_Data.lat_range{i},:,:),[1,2,3,4]);
                moist_means(i)=nanmean(AM4_Data.moist_egr{run_selection}(:,AM4_Data.lat_range{i},:,:),[1,2,3,4]);
            end
            AM4_Data.dry_egr_mean{run_selection} = dry_means;
            AM4_Data.moist_egr_mean{run_selection} = moist_means;
            clear means moist_means dry_means
        end
    end
    save(['AM4_Data.mat'],'AM4_Data')