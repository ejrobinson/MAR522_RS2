%% Batch load CTD V0 is backscatter
tmp_dir = dir('.\CTD'); % IF NEEDED CHANGE DIRECTORY NAME HERE
if exist('t_fold','var')
clear t_fold t_name f_name
end
for i = 1:length(cellfun('isempty',{tmp_dir.name}))
t_fold(i) = string(tmp_dir(i).folder);
t_name(i) = string(tmp_dir(i).name);
end
for i = 1:length(t_name)
    tmp = char(t_name(i));
    if length(tmp) < 3
        df(i) = 1;
    elseif tmp(end-2:end) == 'asc'
        df(i) = 0;
    else
        df(i) = 1;
    end
end
tmp = '\';

t_fold = t_fold(df==0);
t_name = t_name(find(df==0));
slash = repmat(tmp,length(t_name),1);
slash = string(slash);
slash = slash';
f_name = strcat(t_fold,slash,t_name);
clear tmp tmp_dir i df t_fold

CTD_fname = f_name';
for i = 1:length(CTD_fname)
    tmp = char(CTD_fname(i));
    fnames(i) = strcat('ctd_',string(tmp(end-10:end-8)));
    evl_0 = sprintf('%s%s%s%s%s',fnames(i),'=ctd_import(','''',CTD_fname(i),'''',');');
    eval(evl_0)
% ctd_s03 = ctd_import('mar518_s03_RFP.asc');
% ctd_s10 = ctd_import('mar518_s10_RFP.asc');
% ctd_s11 = ctd_import('mar518_s11_RFP.asc');
% ctd_s14 = ctd_import('mar518_s14_RFP.asc');
% ctd_s15 = ctd_import('mar518_s15_RFP.asc');
% ctd_s16 = ctd_import('mar518_s16_RFP.asc');
% ctd_s17 = ctd_import('mar518_s17_RFP.asc');
end
CTD_fname = fnames;
%% Cut CTD to usefull part of profile
%fnames = {'ctd_s03','ctd_s10','ctd_s11','ctd_s14','ctd_s15','ctd_s16','ctd_s17'};
for i = 1:length(fnames)
    evl_0 = sprintf('%s %s','tmp_file =',string(cell2mat(fnames(i))));
    eval(evl_0);
    idx_start = find(tmp_file.s >= 30 & tmp_file.depth > 0.5,1);
    if i == NaN
        idx_end = 620;      
    else
    [~, idx_end] = max(tmp_file.depth);
    end
    vert_vel = idx_start;
    while (tmp_file.depth(vert_vel))-(tmp_file.depth(vert_vel+2)) >= -0.11 & vert_vel < idx_end-2
        vert_vel = vert_vel+1;
    end
    idx_start = vert_vel;
    tmp_names = fieldnames(tmp_file);
    for j = 1:numel(fieldnames(tmp_file))
    evl_1 = sprintf('%s%s%s%s%s%s','tmp_file.',cell2mat(tmp_names(j)),'=','tmp_file.',cell2mat(tmp_names(j)),'(idx_start:idx_end)');
    eval(evl_1)
    end
    evl_2 = sprintf('%s%s',string(cell2mat(fnames(i))), '=tmp_file');
    eval(evl_2)
end
%Housekeeping

clear ans evl_0 evl_1 evl_2 fnames i idx_end idx_start j tmp_file tmp_names vert_vel varname
%% Plot CTD Data
figure(1)
subplot(1,2,1)
hold on
for i = 1:length(CTD_fname) % For each station
    evl_0 = sprintf('%s%s%s%s','plot(',CTD_fname(i),'.s,',CTD_fname(i),'.depth)');
    eval(evl_0);
end
legend(CTD_fname,'location','southeast')
set(gca,'ydir','reverse')
grid on
subplot(1,2,2)
hold on
for i = 1:length(CTD_fname) % For each station
    evl_0 = sprintf('%s%s%s%s','plot(',CTD_fname(i),'.t,',CTD_fname(i),'.depth)');
    eval(evl_0);
end
grid on
legend(CTD_fname,'location','southeast')
set(gca,'ydir','reverse')

%% TriOS data - MUST BE IN SUBFOLDER names TriOS (case sensitve)
root = cd;
% Folder list for TriOS data (f_name)
tmp_dir = dir('.\TriOS'); % IF NEEDED CHANGE DIRECTORY NAME HERE
if exist('t_fold','var')
clear t_fold t_name f_name
end
for i = 1:length(cellfun('isempty',{tmp_dir.name}))
t_fold(i) = string(tmp_dir(i).folder);
t_name(i) = string(tmp_dir(i).name);
end
for i = 1:length(t_name)
    tmp = char(t_name(i));
    if isletter(tmp(1))
        df(i) = 0;
    else
        df(i) = 1;
    end
end
tmp = '\';

t_fold = t_fold(df==0);
t_name = t_name(find(df==0));
slash = repmat(tmp,length(t_name),1);
slash = string(slash);
slash = slash';
f_name = strcat(t_fold,slash,t_name);
clear tmp tmp_dir i df t_fold

% Run though folders and produce Rrs values for each station
for i = 1:length(f_name)
    cur_dir = f_name(i);
    cur_files = dir(sprintf('%s%s',cur_dir,'\*.dat'));
    lowVal = 100;
    
    lowFile = ' ';
    c = 0;
    cd(cur_dir)
    for cFile = 1:length(cur_files)
        fnow = cur_files(cFile).name;
        fid = fopen(fnow,'r');
        for j = 1:43
            fgets(fid);
        end
        for j = 1:196
            tmp_a = fgets(fid);
            tmp_b = textscan(tmp_a,'%f%f%d%d');
            wv(j) = cell2mat(tmp_b(1));
            radiance(j) = cell2mat(tmp_b(2));
        end
        fclose(fid);
        if radiance(50) < lowVal
            figure(2);
            clf;
            plot(wv, radiance)
            %text(320,.2,'Downwelling irradiance : expect 10 or more W/m2/nm; Upwelling radiance : expect 10 or less W/m2/nm');
            %s = input('Enter 1 if good upwelling','s');
            disp('WARNING USER VALIDATION OVERWRITTEN')
            s = '1';
            if s == '1'
                c = cFile;
                lowVal = radiance(50);
                bestRadiance = radiance;
            end
        end
    end
    fnow  = cur_files(c-1).name;
    fid = fopen(fnow,'r');
    for j = 1:43
        fgets(fid);
    end
    for j = 1:196
        tmp_a = fgets(fid);
        tmp_b = textscan(tmp_a,'%f%f%d%d');
        wv(j) = cell2mat(tmp_b(1));
        nearestEdBefore(j) = cell2mat(tmp_b(2));
    end
    fclose(fid);
    fnow = cur_files(c+1).name;
    fid = fopen(fnow,'r');
    for j = 1:43
        fgets(fid);
    end
    for j = 1:196
        tmp_a = fgets(fid);
        tmp_b = textscan(tmp_a,'%f%f%d%d');
        wv(j) = cell2mat(tmp_b(1));
        nearestEdAfter(j) = cell2mat(tmp_b(2));
    end
    fclose(fid);
%     figure(2); clf;
%     plot(wv, nearestEdBefore, 'k','linewidth',2);
%     hold on
%     plot(wv, nearestEdAfter, 'k--');
%     plot(wv, bestRadiance, 'g','linewidth',2);
%     set(gca,'yscale','log');
%     legend({'Ed before', 'Ed after', 'lowest Lu'});
%     grid on;
    Rrs_1 = bestRadiance ./ (nearestEdBefore.*(2*pi));    % THESE MAY NEED TO CHANGE 
    Rrs_2 = bestRadiance ./ (nearestEdAfter.*(2*pi));     % SEE NOTE P4 Practical guide 
% 
%     figure(2); clf;
%     plot(wv, Rrs_1, 'k','linewidth',2)
%     hold on
%     plot(wv, Rrs_2, 'k--');
    evl_0 = sprintf('%s%s%s','Rrs_',t_name(i),'=Rrs_2;');
    eval(evl_0);
    evl_0 = sprintf('%s%s%s%s%s','Rrs_fname(i) =string(','''', 'Rrs_',t_name(i),'''',');');
    eval(evl_0);
    evl_1 = sprintf('%s%s%s%s%s%s','save(','''','Rrs_',t_name(i),'''',',','''','Rrs_',t_name(i),'''',');');
    eval(evl_1);
end
cd(root)
%Housekeeping

clear ans bestRadiance c cFile cur_dir cur_files evl_0 evl_1 f_name fid fnow i j lowFile lowVal nearestEdAfter nearestEdBefore radiance Rrs_1 Rrs_2 s slash t_name tmp_a tmp_b 

%% Plot TriOS data on one plot
figure(3)
clf
%hold on
%plot(wv,Rrs_stn03);
%plot(wv,Rrs_stn10);
%plot(wv,Rrs_stn11);
%plot(wv,Rrs_stn14);
%plot(wv,Rrs_stn15);
%plot(wv,Rrs_stn16);
%plot(wv,Rrs_stn17);
%legend(CTD_fname,'location','northeast')
hold on
for i = 1:length(Rrs_fname) % For each station
    evl_0 = sprintf('%s%s%s','plot(wv,',Rrs_fname(i),')');
    eval(evl_0);
end
grid on
legend(CTD_fname,'location','southeast')

%% TriOS Band exctraction
% Sentinel 2 red band 649.1-680.1
rb4s = 541.8;
rb4e = 577.8;
b4_idx_start = find(wv >= rb4s,1);
b4_idx_end = find(wv >= rb4e,1) -1;
% Exctract red band means
for i = 1:length(Rrs_fname)
    evl_0 = sprintf('%s%s','tmp =',Rrs_fname(i)) 
    eval(evl_0)
    tmp = tmp(b4_idx_start:b4_idx_end);
    tmp_mean_b4 = nanmean(tmp);
    evl_1 = sprintf('%s%s',Rrs_fname(i),'_b4m=tmp_mean_b4;')
    eval(evl_1);
    evl_1 = sprintf('%s%s%s%s%s','Rrs_b4_name(i) =string(','''',Rrs_fname(i),'_b4m','''',');');
    eval(evl_1);
end
% Sentinel 2 red band 649.1-680.1
nir_s = 779.8;
nir_e = 885.8
nir_idx_start = find(wv >= nir_s,1);
nir_idx_end = find(wv >= nir_e,1) -1;
% Exctract red band means
for i = 1:length(Rrs_fname)
    evl_0 = sprintf('%s%s','tmp =',Rrs_fname(i)) 
    eval(evl_0)
    tmp = tmp(nir_idx_start:nir_idx_end);
    tmp_mean_nir = nanmean(tmp);
    evl_1 = sprintf('%s%s',Rrs_fname(i),'_nir=tmp_mean_nir;')
    eval(evl_1);
    evl_1 = sprintf('%s%s%s%s%s','Rrs_nir_name(i) =string(','''',Rrs_fname(i),'_nir','''',');');
    eval(evl_1);
end
%Housekeeping

clear evl_0 evl_1 tmp_mean_b4 b4_idx_start b4_idx_end
%% PRR File Import
tmp_dir = dir('.\PRR'); % IF NEEDED CHANGE DIRECTORY NAME HERE
if exist('t_fold','var')
clear t_fold t_name f_name
end
for i = 1:length(cellfun('isempty',{tmp_dir.name}))
t_fold(i) = string(tmp_dir(i).folder);
t_name(i) = string(tmp_dir(i).name);
end
for i = 1:length(t_name)
    tmp = char(t_name(i));
    if length(tmp) < 3
        df(i) = 1;
    elseif tmp(end-2:end) == 'DAT'
        df(i) = 0;
    else
        df(i) = 1;
    end
end
tmp = '\';

t_fold = t_fold(df==0);
t_name = t_name(find(df==0));
slash = repmat(tmp,length(t_name),1);
slash = string(slash);
slash = slash';
f_name = strcat(t_fold,slash,t_name);
clear tmp tmp_dir i df t_fold

PRR_fname = f_name';
for i = 1:length(PRR_fname)
    tmp = char(PRR_fname(i));
    fnames(i) = lower(strcat('prr_',string(tmp(end-6:end-4))));
    evl_0 = sprintf('%s%s%s%s%s',fnames(i),'=prr_import(','''',PRR_fname(i),'''',');');
    eval(evl_0);
end
PRR_fname = fnames;
%Housekeeping

clear evl_0 f_name fnames i slash t_name tmp

%% PRR Cleaning

%For each PRR Station
for i = 1:length(PRR_fname)
    % Copy current station to tmp_wd
    evl_0 = sprintf('%s%s','tmp_wd=',PRR_fname(i));
    eval(evl_0)
    id_str = fieldnames(tmp_wd);
    % For each vector list the names
    for j = 1:length(id_str)
        tmp = char(id_str(j));
        if length(tmp) < 2
            df(j) = 1;
        elseif tmp(1:2) == 'Ed'
            df(j) = 0;
        else
            df(j) = 1;
        end
    end
    clear tmp
    % Extract those begining with Ed (channels) to ch_n
    ch_n = id_str(find(df == 0));
    ch_n = string(ch_n);
    ch_n = ch_n(1:end-2);
    for j = 1:length(ch_n)
        tmp_ch = char(ch_n(j));
    ch_val(j) = str2double(tmp_ch(end-2:end));
    end
    %Housekeeping

    clear tmp_ch
    % Calculate start and finish indexes [idx_start idx_end]
        % Select start value as when depth exceded 0.6m minus 4 indexes based of most samples hovering above 0.6m for stabalisation, 4 point saftey buffer
    tmp_d = tmp_wd.d;
    idx_startp = find(tmp_d >= 0.6,1);
    if idx_startp >4
        idx_start = idx_startp - 4;
    else
        idx_start = idx_startp;
    end
    [~,idx_end] = max(tmp_d);
    tmp_d_c = tmp_d(idx_start:idx_end);
    tmp_wd.d_c = tmp_d_c;
    % For each channel (i_l j_f)
    for j = 1:length(ch_n)
        % Pull to temporary storage tmp_d and tmp_ch (i_l j_l)
        evl_1 = sprintf('%s%s%s','tmp_ch = tmp_wd.',ch_n(j),';');
        eval(evl_1)

        % Cut to start and end index (i_l j_l) 
        tmp_ch = tmp_ch(idx_start:idx_end);
        % Copy back into struct with _c append for clean data(i_l j_l)
        evl_2 = sprintf('%s%s%s','tmp_wd.',ch_n(j),'_c=tmp_ch;');
        eval(evl_2)
        cd_names(j) = string(sprintf('%s%s',ch_n(j),'_cln'));
        negvals = find(tmp_ch <= 0);
        tmp_ch(negvals) = NaN;
        tmp_ch_ln = log(tmp_ch);
        evl_2 = sprintf('%s%s%s','tmp_wd.',ch_n(j),'_cln=tmp_ch_ln;');
        eval(evl_2)
        % Calculare p(m) and r values(i_l j_l)
        % m value is (i_l j_l)
        [tmp_pfit, r_val]= polyfit(tmp_d_c,tmp_ch_ln,1);
        tmp_m(j) = -tmp_pfit(1);
        %tmp_r = r_val.normr;
        tmp_r(j) = 1 - (r_val.normr/norm(tmp_ch_ln - mean(tmp_ch_ln)))^2;
        evl_2 = sprintf('%s%s%s','tmp_wd.',ch_n(j),'_r=tmp_r(j);');
        eval(evl_2)
        evl_2 = sprintf('%s%s%s','tmp_wd.',ch_n(j),'_m=tmp_m(j);');
        eval(evl_2)
    end
    % Copy all r & kD vals into one variable (i_l j_f)
    tmp_wd.r2_vals = tmp_r;
    tmp_wd.kD_vals = tmp_m;
    % Copy Back to base workspace (i_l j_f)
    evl_3 = sprintf('%s%s',PRR_fname(i),'=tmp_wd;');
    eval(evl_3);
end
%Housekeeping
clear df evl_0 evl_1 evl_2 evl_3 i id_str idx_end idx_start idx_startp j negvals tmp_ch tmp_ch_ln tmp_d tmp_d_c tmp_wd tmp_r tmp_m r_val tmp_pfit

%% Plotting all kD vals (i free j free)
figure(4)
hold on
for i = 1:length(PRR_fname) % For each station
    evl_0 = sprintf('%s%s%s','plot(ch_val,',PRR_fname(i),'.kD_vals)');
    eval(evl_0);
end
grid on
set(get(gca, 'YLabel'), 'String', 'kD m^-^1');
set(get(gca, 'XLabel'), 'String', 'PRR Band');
legend(PRR_fname,'location','southeast');
xticks(ch_val)
xticklabels(ch_n)

clear evl_0 i
%% Conversion of Rrs into Kd PRR ch5 =Sentinel band 4 (not currently set to this as changed)
% Cancatonate values so can actually fit the trend line becasue matlab is stupid sometimes
for i = 1:length(PRR_fname)
evl_0 = sprintf('%s%s%s', 'kd_vals_5(i) =',PRR_fname(i),'.kD_vals(2);'); %% PRR BAND SELECT
eval(evl_0 )
evl_1 = sprintf('%s%s%s','rrs_vals_5(i)=',Rrs_b4_name(i),';');
eval(evl_1);
end
figure(5)
subplot(2,2,1)
scatter(rrs_vals_5,kd_vals_5);
tsm_vals = [0.00436 0.00307 0.00175 0.00124 0.001716667 0.00182 0.001896667 0.001456667 0.002146667 0.002753333 0.001133333 0.001173333 0.00164];
subplot(2,2,2)
scatter(rrs_vals_5,tsm_vals);

%% Ratio calc
for i = 1:length(Rrs_b4_name)
    evl_0 = sprintf('%s%s%s%s','rrs_ratio(i) = ',Rrs_nir_name(i),'/',Rrs_b4_name(i))
    eval(evl_0)
end
%% Sort vals
rrs_kd = [rrs_vals_5', kd_vals_5'];
rrs_tsm = [rrs_vals_5', tsm_vals'];
rrs_tsm_ratio = [rrs_ratio' tsm_vals'];


[~,idx] = sort(rrs_kd(:,1));
rrs_kd_s = rrs_kd(idx,:);

[~,idx] = sort(rrs_tsm(:,1));
rrs_tsm_s = rrs_tsm(idx,:);

[~,idx] = sort(rrs_tsm_ratio(:,1));
rrs_tsmr_s = rrs_tsm_ratio(idx,:);

[pfit_rrs_tsm, rval_rrs_tsm] = polyfit(rrs_tsm_s(:,1),rrs_tsm_s(:,2),1);
[pfit_rrsr_tsm, rval_rrsr_tsm] = polyfit(rrs_tsmr_s(:,1),rrs_tsmr_s(:,2),1);
[pfit_rrs_kd, rval_rrs_kd] = polyfit(rrs_kd_s(:,1),rrs_kd_s(:,2),1);

 fit_rrs_tsm = fit(rrs_tsm_s(:,1),rrs_tsm_s(:,2),'poly1');
 fit_rrsr_tsm = fit(rrs_tsmr_s(:,1),rrs_tsmr_s(:,2),'poly1');
 fit_rrs_kd = fit(rrs_kd_s(:,1),rrs_kd_s(:,2),'poly1');
% 
 confint_tsm = confint(fit_rrs_tsm,0.95);
 confint_kd = confint(fit_rrs_kd,0.95);
 confint_rrsr = confint(fit_rrsr_tsm,0.95);

% 
 predb_tsm = predint(fit_rrs_tsm,rrs_tsm_s(:,1));
 predb_kd = predint(fit_rrs_kd,rrs_tsm_s(:,1));
 predb_rrsr = predint(fit_rrsr_tsm,rrs_tsmr_s(:,1));

r_tsm = 1- (rval_rrs_tsm.normr/norm(tsm_vals - mean(tsm_vals)))^2;
r_tsr = 1- (rval_rrsr_tsm.normr/norm(tsm_vals - mean(tsm_vals)))^2;
r_kd  = 1- (rval_rrs_kd.normr/norm(kd_vals_5 - mean(kd_vals_5)))^2;        

pred_tsm = polyval([fit_rrs_tsm.p1 fit_rrs_tsm.p2],rrs_tsm_s(:,1));
pred_rrsr = polyval([fit_rrsr_tsm.p1 fit_rrsr_tsm.p2],rrs_tsmr_s(:,1));
pred_kd = polyval([fit_rrs_kd.p1 fit_rrs_kd.p2],rrs_tsm_s(:,1));



 figure(5)
 subplot(2,2,2)
 hold on
 plot(rrs_tsm_s(:,1),predb_tsm,'r--')
 plot(rrs_tsm_s(:,1),pred_tsm,'r')
set(get(gca, 'XLabel'), 'String', 'Mean Rrs (541.8-577.8um) mW/(m^2 nm Sr)')
set(get(gca, 'YLabel'), 'String', 'Total suspended matter (g)');
grid on

 subplot(2,2,1)
 hold on
 plot(rrs_tsm_s(:,1),predb_kd,'r--')
 plot(rrs_tsm_s(:,1),pred_kd,'r-')
 set(get(gca, 'XLabel'), 'String', 'Mean Rrs (541.8-577.8um) mW/(m^2 nm Sr)')
 set(get(gca, 'YLabel'), 'String', 'Kd (m^-^1)');
 clear evl_0 evl_1 r4be r4bs i
 grid on
 
 subplot(2,2,4)
 hold on
 scatter(kd_vals_5,tsm_vals)
 set(get(gca, 'XLabel'), 'String', 'Kd (m^-^1)');
 set(get(gca, 'YLabel'), 'String', 'Total suspended matter (g)');
grid on

%% Compute Ratio of RRS(b4) and RRS(nir)

subplot(2,2,3)
hold on
scatter(rrs_ratio,tsm_vals)
 plot(rrs_tsmr_s(:,1),predb_rrsr,'r--')
 plot(rrs_tsmr_s(:,1),pred_rrsr,'r')
set(get(gca, 'XLabel'), 'String', 'Rrs ratio (rrs(833)/rrs(560)') %% NEED TO EXPLAIN THESE NOT EXACT 
set(get(gca, 'YLabel'), 'String', 'Total suspended matter (g)');