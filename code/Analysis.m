%%%%%% Analysis of airfoil (A) - flag (F) setup with %%%%%%%
%%%%%% hydrodynamic and electromechanical coupling   %%%%%%%

format long
clear all; close all; clc;

FTsz = 20; 
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

%% Parameters
flag_save = -1; % save LDV downsampled ts
flag_norm = 1; % >0: normalize TE and cond TE, <0: no normalization
frame_rate=60; % No. of frames per second captured by camera - for motion of A & F
down_fps=30; % No. of frames per second timeseries should be downsampled
             % frame_skip (delf) = frame_rate/down_fps
flag_sa = 1; % >0 : seasonally adjust data, <0: do not adjust data
upsample=4; % upsampling factor for resampling LDV timeseries to match 
            % in time with A & F

m_embed = 2; % embedding dimension for symbolizing the time series
Ndel = down_fps; % max delay upto a second
delarr = 1:1:Ndel; 
Nsurr = 10000; % No. of samples for surrogate distribution

tidx = 4401:5100; % portion of the time series to plot

%% Data
data_dir = './../data/';
str_delay = '03'; % delay of mech coupling of flag wrt a/f: 00,01,02,03
str_chanfreq = '19hz'; % water channelf requency: 0hz, 3hz, 19hz
str_ldvloc = 'ct'; % location in channel where ldv measurement is taken: ct
str_flagtype = 's2'; % type of flag --> s1: passive flag, s2: active flag 
str_noise = 'no'; % 'no':noise, '':no noise
str_Delt_noise = '_d008'; % mean time interval between noise startles in sec:
                         % '': if no noise, '_d004':0.04s, '_d008':0.08s
filename = [str_delay '_' str_chanfreq '_' str_ldvloc '_' ...
            str_flagtype str_noise str_Delt_noise];


%% Other directories
files_dir = data_dir;
save_fig_dir = './figures/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end

%% Load files
% Read A and F data
load([files_dir, filename, '_Flag2.mat']); % 'time', 'flagtipY'
load([files_dir, filename, '_Foil.mat']); % 'time', 'foilang'
% Read LDV file (both raw and already moving averaged data)
load([files_dir, filename, '_LDV.mat']); %, 'tmov', 'umov', 'tLDV', 'uLDV' % we are not using the moving avg u

flagfrontY_raw = flagfrontY;
flagtipY_raw = flagtipY;
foilang_raw = foilang;

%% Downsample/resample LDV's nonuniform timeseries 
 % to a uniform one concurrent with A & F time series

fs = down_fps; %1/ts; % uniform sampling rate
[umov,tmov] = resample(uLDV,tLDV,fs,upsample,1); % upsample and then downsample to fs

if (flag_save>0)
    if ~exist(files_dir,'dir')
        mkdir(files_dir);
    end
    savefname = [files_dir, filename, '_LDVdowns.mat'];
    save(savefname, 'tmov', 'umov');
end

%% Seasonally adjust
if (flag_sa>0)
    disp('Read seasonally adjusted time series of A and F ...');
    % read from python MSTL saved files
    load([files_dir, filename, '_Flag2_SA.mat']); % 'flagtipY_trend','flagtipY_seas','flagtipY_resid'
    load([files_dir, filename, '_Foil_SA.mat']); % 'foilang_trend','foilang_seas','foilang_resid'
    flagfrontY = (flagfrontY_trend + flagfrontY_resid)';
    flagtipY = (flagtipY_trend + flagtipY_resid)';
    foilang = (foilang_trend + foilang_resid)';
end

%% downsample to -- A (airfoil), F (flag)
time = time'; 
[F2,tarr] = downsample(flagfrontY,time,down_fps,frame_rate); 
[F,tarr] = downsample(flagtipY,time,down_fps,frame_rate); 
[A,tarr] = downsample(foilang,time,down_fps,frame_rate); 

% Match the times
if (size(tarr,1)<size(tmov,1)) 
    tmov=tmov(1:size(tarr,1)); 
    umov=umov(1:size(tarr,1));
elseif (size(tarr,1)>size(tmov,1)) 
    tarr=tarr(1:size(tmov,1)); 
    A=A(1:size(tmov,1));
    F=F(1:size(tmov,1));
    F2=F2(1:size(tmov,1));
end
if (abs(tarr-tmov)>1e-4)
    disp('Error! LDV and flag/foil ts do not match in time after downsampling!');
    disp(tarr(1:5));disp(tmov(1:5));
    pause;
else
    t = tarr;
end

u = umov;

%% Plot all the processed time series

dt = (t(2)-t(1));
tst = ceil(tidx(1)*dt/5)*5; % round to nearest higher 5
tend = floor(tidx(end)*dt/5)*5; % round to nearest lower 5

f1=figure;
subplot(4,1,1);
% Airfoil
plot(t(tidx),A(tidx), '.','Color',[0 0.4470 0.7410],'MarkerSize',12); hold on; %title('Airfoil angle (deg.)');
% Plot raw ts for same time duration to compare
st_time = min(t(tidx)); end_time = max(t(tidx));
[st_time_A,st_idx_A] = min(abs(time-st_time));
[end_time_A,end_idx_A] = min(abs(time-end_time));
plot(time(st_idx_A:end_idx_A), foilang(st_idx_A:end_idx_A), ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([-50 50]);
ylabel('A');
% Flag
subplot(4,1,2); 
plot(t(tidx),F2(tidx)*100, '.','Color',[0.8500 0.3250 0.0980],'MarkerSize',12); hold on; %title('Flag front deflection (cm)'); 
% Plot raw ts for same time duration to compare
plot(time(st_idx_A:end_idx_A), flagfrontY(st_idx_A:end_idx_A)*100, ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([-0.6 0.6]); yticks([-0.6 0 0.6]);
ylabel('F_L'); % flag front
% Flag2
subplot(4,1,3); 
plot(t(tidx),F(tidx)*100, '.','Color',[0.4660 0.6740 0.1880],'MarkerSize',12); hold on; %title('Flag deflection (cm)'); 
% Plot raw ts for same time duration to compare
plot(time(st_idx_A:end_idx_A), flagtipY(st_idx_A:end_idx_A)*100, ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([-5 5]); yticks([-5 0 5]);
ylabel('F_T'); % flag tip
% LDV
subplot(4,1,4); 
plot(t(tidx),u(tidx), '.','Color',[0.3010 0.7450 0.9330],'MarkerSize',12); hold on; %title('LDV velocity (m/s)'); 
% Plot raw ts of LDV for same time duration to compare
st_time = min(t(tidx)); end_time = max(t(tidx));
[st_time_LDV,st_idx_LDV] = min(abs(tmov-st_time));
[end_time_LDV,end_idx_LDV] = min(abs(tmov-end_time));
plot(tmov(st_idx_LDV:end_idx_LDV), umov(st_idx_LDV:end_idx_LDV), ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([0.1 0.7]); yticks([0.1 0.7]);
ylabel('u');
xlabel('time (s)');

f1.Position = [1 100 564 540];
saveas(gcf,[save_fig_dir filename '_processed_ts' '_ups' num2str(upsample)],'png');


%% Plot Raw time series
f0=figure;
subplot(4,1,1);
% Airfoil
st_time = min(t(tidx)); end_time = max(t(tidx));
[~,st_idx_A] = min(abs(time-st_time));
[~,end_idx_A] = min(abs(time-end_time));
plot(time(st_idx_A:end_idx_A), foilang_raw(st_idx_A:end_idx_A), ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([-50 50]);
ylabel('\theta (deg.)');
% Flag
subplot(4,1,2);
plot(time(st_idx_A:end_idx_A), flagfrontY_raw(st_idx_A:end_idx_A)*100, ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([-0.6 0.6]); yticks([-0.6 0 0.6]);
ylabel('\Delta y_L (cm)');
% Flag
subplot(4,1,3);
plot(time(st_idx_A:end_idx_A), flagtipY_raw(st_idx_A:end_idx_A)*100, ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([-5 5]); yticks([-5 0 5]);
ylabel('\Delta y_T (cm)');
% LDV
subplot(4,1,4); 
st_time = min(t(tidx)); end_time = max(t(tidx));
[~,st_idx_LDV] = min(abs(tLDV-st_time));
[~,end_idx_LDV] = min(abs(tLDV-end_time));
% plot(t(tidx),u(tidx), '-om','LineWidth',1.5,'MarkerSize',5); hold on; %title('LDV velocity (m/s)'); 
plot(tLDV(st_idx_LDV:end_idx_LDV), uLDV(st_idx_LDV:end_idx_LDV), ...
          '-k','LineWidth',1.5);
xlim([tst tend]); ylim([0.1 0.7]); yticks([0.1 0.7]);
ylabel('U (m/s)');
xlabel('time (s)');

f0.Position = [1 100 564 540];
saveas(gcf,[save_fig_dir filename '_raw_ts' '_ups' num2str(upsample)],'png');


%% Symbolize time series
piA = symbolize_data(A,m_embed);
piF = symbolize_data(F,m_embed);
piF2 = symbolize_data(F2,m_embed);
piu = symbolize_data(u,m_embed);


%% Plot symbolized time series

f8=figure;
subplot(4,1,1);
% Airfoil
plot(t(tidx),piA(tidx), 'Color',[0 0.4470 0.7410],'LineWidth',1.5,'MarkerSize',12); hold on; %title('Airfoil angle (deg.)');
ylabel('\pi_A'); xlim([tst tend]); 
% Flag
subplot(4,1,2); 
plot(t(tidx),piF2(tidx), 'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5,'MarkerSize',12); hold on; %title('Flag front deflection (cm)'); 
ylabel('\pi_{F_L}'); xlim([tst tend]); 
% Flag2
subplot(4,1,3); 
plot(t(tidx),piF(tidx), 'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'MarkerSize',12); hold on; %title('Flag tip deflection (cm)'); 
ylabel('\pi_{F_T}'); xlim([tst tend]); 
% LDV
subplot(4,1,4); 
plot(t(tidx),piu(tidx), 'Color',[0.3010 0.7450 0.9330],'LineWidth',1.5,'MarkerSize',12); hold on; %title('LDV velocity (m/s)'); 
ylabel('\pi_u'); xlim([tst tend]); 
xlabel('time (s)');

f8.Position = [1 100 564 540];
saveas(gcf,[save_fig_dir filename '_symb_ts' '_ups' num2str(upsample)],'png');


%% Transfer entropy

% Calc TE_A->F and TE_F->A
TA_F = nan(Ndel,1); TF_A = nan(Ndel,1);
for idel = 1:Ndel
    del = delarr(idel);
    TA_F(idel) = transfer_entropy_delay(piF,piA,del,flag_norm);
    TF_A(idel) = transfer_entropy_delay(piA,piF,del,flag_norm);    
end
netTA_F = TA_F - TF_A;

% Plot 
f2=figure;
plot(delarr,TA_F,'-r','LineWidth',2); hold on;
plot(delarr,TF_A,'-b','LineWidth',2); 
legend('TE_{A \rightarrow {F_T}}','TE_{{F_T} \rightarrow A}');
xlabel('\delta');
ylim([0 0.05]); yticks([0.0 0.01 0.02 0.03 0.04 0.05]);
f2.Position = [750 400 500 250];
saveas(gcf,[save_fig_dir filename '_delayTE' '_ups' num2str(upsample)],'png');

%% Pick the delay with max TE_A_to_F 
[maxTA_F,maxTA_F_idx] = max(TA_F);
del = delarr(maxTA_F_idx)
[maxnetTA_F,maxnetTA_F_idx] = max(netTA_F);

TA_Fval = maxTA_F
TF_Aval = TF_A(maxTA_F_idx);

%% Shuffle source time series and test statistical significance
TA_Fsurr = zeros(Nsurr,1);
for i=1:Nsurr
    [piAshuff,~,~] = permutate(piA,piF,ones(size(piA)));
    TA_Fsurr(i) = transfer_entropy_delay(piF,piAshuff,del,flag_norm);
end
TA_Fcutoff = prctile(TA_Fsurr,95);

TA_Fsurr = [TA_Fsurr; TA_Fval];
TA_Fsurr = sort(TA_Fsurr);
pval = 1-(find(TA_Fsurr==TA_Fval)/length(TA_Fsurr))

f4=figure; 
h=histogram(TA_Fsurr,'Normalization','pdf');
edges = h.BinEdges; 
pdf = h.Values;
centers = 0.5*(edges(1:end-1)+edges(2:end));
semilogx(centers,pdf,'-k','LineWidth',2); hold on;
semilogx(TA_Fcutoff*ones(100,1),linspace(0,max(pdf)+500,100),'--r','LineWidth',2); hold on;
semilogx(TA_Fval*ones(100,1),linspace(0,max(pdf)+500,100),'-g','LineWidth',2); 
xlabel('TE_{A \rightarrow {F_T}}'); ylabel('PDF'); 
f4.Position = [350 400 400 250];
xlim([1e-5 1e-1]); xticks([1e-5 1e-4 1e-3 1e-2 1e-1]);
saveas(gcf,[save_fig_dir filename '_TE_AF' '_ups' num2str(upsample)],'png');


%% Pick the delay with max TE_u_to_F 

% Calc TE A->F and F->A symbolic
Tu_F = nan(Ndel,1); 
for idel = 1:Ndel
    del_tmp = delarr(idel);
    Tu_F(idel) = transfer_entropy_delay(piF,piu,del_tmp,flag_norm);
end
[maxTu_F,maxTu_F_idx] = max(Tu_F);
deluF = delarr(maxTu_F_idx)
Tu_Fval = maxTu_F;

%% Cond TE_ A->F | u

TA_F_u = nan(Ndel,1);
for idel = 1:Ndel
    del2 = delarr(idel);
    TA_F_u(idel) = cond_transfer_entropy_delay(piF,piA,piu,del,del2,flag_norm); 
end
    
TA_F_u_deluF = cond_transfer_entropy_delay(piF,piA,piu,del,deluF,flag_norm); 

% Plot 
f6=figure;
plot(delarr(1:del),TA_Fval*ones(size(delarr(1:del),2),1),'--k','LineWidth',2); hold on;
plot(delarr(1:del),TA_F_u(1:del),'-r','LineWidth',2); hold on; 
xlabel('\delta_2');
ylabel('TE_{A \rightarrow {F_T}|u}');
f6.Position = [950 100 320 250];
saveas(gcf,[save_fig_dir filename '_delaycondTE_AFu' '_ups' num2str(upsample)],'png');

% Pick the delay with min cond TE 
[minTA_F_u,minTA_F_u_idx] = min(TA_F_u(2:del-1));
del2 = delarr(minTA_F_u_idx)
TA_F_uval = minTA_F_u

% Shuffle u (cond variable) maintaining the A-F dynamics to test significance  
TA_F_usurr = zeros(Nsurr,1);
for i=1:Nsurr
    [piushuff,~,~] = permutate(piu,piA,piF);
    TA_F_usurr(i) = cond_transfer_entropy_delay(piF,piA,piushuff,del,del2,flag_norm);
end
TA_F_ucutoff = prctile(TA_F_usurr,5);

TA_F_usurr = [TA_F_usurr; TA_F_uval];
TA_F_usurr = sort(TA_F_usurr);
pval = (find(TA_F_usurr==TA_F_uval)/length(TA_F_usurr))

f7=figure; 
h=histogram(TA_F_usurr,'Normalization','pdf');
edges = h.BinEdges; 
pdf = h.Values;
centers = 0.5*(edges(1:end-1)+edges(2:end));
semilogx(centers,pdf,'-k','LineWidth',2); hold on;
semilogx(TA_F_ucutoff*ones(100,1),linspace(0,max(pdf)+500,100),'--r','LineWidth',2); hold on;
semilogx(TA_F_uval*ones(100,1),linspace(0,max(pdf)+500,100),'-g','LineWidth',2); 
xlabel('TE_{A \rightarrow {F_T}|u}'); ylabel('PDF'); 
f7.Position = [350 400 400 250];
saveas(gcf,[save_fig_dir filename '_condTE_AFu' '_ups' num2str(upsample)],'png');


%% Pick the delay with max TE_F2_to_F 
TF2_F = nan(Ndel,1); 
for idel = 1:Ndel
    del_tmp = delarr(idel);
    TF2_F(idel) = transfer_entropy_delay(piF,piF2,del_tmp,flag_norm);
end
[maxTF2_F,maxTF2_F_idx] = max(TF2_F);
delF2F = delarr(maxTF2_F_idx)
TF2_Fval = maxTF2_F;

%% Cond TE_ A->F | F2

TA_F_F2 = nan(Ndel,1);
for idel = 1:Ndel
    del2 = delarr(idel);
    TA_F_F2(idel) = cond_transfer_entropy_delay(piF,piA,piF2,del,del2,flag_norm); 
end
TA_F_F2_delF2F = cond_transfer_entropy_delay(piF,piA,piF2,del,delF2F,flag_norm); 

% Plot 
f6=figure;
plot(delarr(1:del),TA_Fval*ones(size(delarr(1:del),2),1),'--k','LineWidth',2); hold on;
plot(delarr(1:del),TA_F_F2(1:del),'-r','LineWidth',2); 
xlabel('\delta_2'); 
ylabel('TE_{A \rightarrow {F_T}|{F_L}}');
f6.Position = [950 100 320 250];
saveas(gcf,[save_fig_dir filename '_delaycondTE_AF_F2' '_ups' num2str(upsample)],'png');

% Pick the delay with min cond TE 
[minTA_F_F2,minTA_F_F2_idx] = min(TA_F_F2(1:del-1)); 
del2 = delarr(minTA_F_F2_idx)
TA_F_F2val = minTA_F_F2;

% Shuffle u (conditioning variable) maintaining the A-F dynamics to test significance  
TA_F_F2surr = zeros(Nsurr,1);
for i=1:Nsurr
    [piF2shuff,~,~] = permutate(piF2,piA,piF);
    TA_F_F2surr(i) = cond_transfer_entropy_delay(piF,piA,piF2shuff,del,del2,flag_norm);
end
TA_F_F2cutoff = prctile(TA_F_F2surr,5);

TA_F_F2surr = [TA_F_F2surr; TA_F_F2val];
TA_F_F2surr = sort(TA_F_F2surr);
pval = (find(TA_F_F2surr==TA_F_F2val)/length(TA_F_F2surr))

f7=figure; 
h=histogram(TA_F_F2surr,'Normalization','pdf');
edges = h.BinEdges; 
pdf = h.Values;
centers = 0.5*(edges(1:end-1)+edges(2:end));
semilogx(centers,pdf,'-k','LineWidth',2); hold on;
semilogx(TA_F_F2cutoff*ones(100,1),linspace(0,max(pdf)+500,100),'--r','LineWidth',2); hold on;
semilogx(TA_F_F2val*ones(100,1),linspace(0,max(pdf)+500,100),'-g','LineWidth',2); 
xlabel('TE_{A \rightarrow {F_T}|{F_L}}'); ylabel('PDF'); 
f7.Position = [350 400 400 250];
saveas(gcf,[save_fig_dir filename '_condTE_AF_F2' '_ups' num2str(upsample)],'png');