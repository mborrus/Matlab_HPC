% Data Job

clear
cd '/home/users/mborrus/Matlab_HPC'
addpath('/home/users/mborrus/Matlab_HPC/scripts/EGR')
addpath('/home/users/mborrus/Matlab_HPC/plots')
addpath('/home/users/mborrus/Matlab_HPC/scripts/cbrewer')
load('./data/EGR/lambda.mat')
load('AM4_Data.mat')

for run_selection = [4:9];
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
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff); clear var_RMS var_Mean var_STD
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp);
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
    for run_selection = [4:9];
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
            var_long(i,:,:,:,:) = temp_var(:,:,1:24,1:100);
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Np,Ntime] = size(var_diff);
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp);
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
    for run_selection = [4:9];
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
            var_long(i,:,:,:) = temp_var(:,:,1:100);
            end
        end
        clear var_1 temp_var
        
        [Nrun,Nlon,Nlat,Ntime] = size(var_diff);
        for lats = 1:6
            for pres=1:3
                range_temp = AM4_Data.lat_range{lats}; rng = length(range_temp);
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