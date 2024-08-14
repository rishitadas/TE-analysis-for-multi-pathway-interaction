%%%%%% Analysis of airfoil (A) - flag (F) setup with %%%%%%%
%%%%%% hydrodynamic coupling only                    %%%%%%%

format long
clear all;
clc;
close all;
FTsz = 20; 
set(groot,'defaultAxesFontSize',FTsz);
set(groot,'defaultLegendFontSize',FTsz);

%% Parameters
flag_save = 1; % save LDV downsampled ts
flag_norm = 1; % >0: normalize TE and cond TE, <0: no normalization
frame_rate=60; % No. of frames per second captured by camera - for motion of A & F
down_fps=30; % No. of frames per second timeseries should be downsampled
             % frame_skip (delf) = frame_rate/down_fps
%tavg = 0.005; tavg_str='0_005'; % avg over this time distance -- keep a bit smaller than delt
                 % e.g. 0.035 for movavg_fps=30 was used before
flag_sa = 1; % >0 : seasonally adjust data, <0: do not adjust
upsample=4; % upsampling factor for resampling LDV timeseries to match 
            % in time with A & F

m_embed = 2; % embedding dimension for symbolization
Ndel = down_fps; % Max delay upto a second
delarr = 1:1:Ndel; 
Nsurr = 10000; % No. of samples for surrogate distribution

tidx = 3401:3900; % To plot time series
ttot = tidx/down_fps;
tminplot = ceil(ttot(1)/5)*5;
tmaxplot = floor(ttot(end)/5)*5;

%% Data  
data_dir = './../data/';
str_delay = '00'; % delay of mech coupling of flag wrt a/f: 00,01,02,03
str_chanfreq = '19hz'; % water channelf requency: 0hz, 3hz, 19hz
str_ldvloc = 'ct'; % location in channel where ldv measurement is taken: ct
str_flagtype = 's1'; % type of flag --> s1: passive flag, s2: active flag 
str_noise = ''; % 'no':noise, '':no noise
str_Delt_noise = ''; % mean time interval between noise startles in sec:
                         % '': if no noise, '_d004':0.04s, '_d008':0.08s
filename = [str_delay '_' str_chanfreq '_' str_ldvloc '_' ...
            str_flagtype str_noise str_Delt_noise];

%% other directories
files_dir = data_dir;
save_fig_dir = './figures/';
if exist(save_fig_dir, 'dir')==0
    mkdir(save_fig_dir);
end

%% Load files
% Read flag and foil data
load([files_dir, filename, '_Flag2.mat']); % 'time', 'flagtipY'
load([files_dir, filename, '_Foil.mat']); % 'time', 'foilang'
% Read LDV file (both raw and already moving averaged data)
load([files_dir, filename, '_LDV.mat']); %, 'tmov', 'umov', 'tLDV', 'uLDV' % we are not using the moving avg u

flagtipY_raw = flagtipY;
foilang_raw = foilang;

%% Downsample/resample LDV's nonuniform timeseries 
% to a uniform one concurrent with A & F time series
fs = down_fps; %1/ts; % uniform sampling rate
[umov,tmov] = resample(uLDV,tLDV,fs,upsample,1); % upsample by 3 and then downsample to fs

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
    load([files_dir, filename, '_Flag_SA.mat']); % 'flagtipY_trend','flagtipY_seas','flagtipY_resid'
    load([files_dir, filename, '_Foil_SA.mat']); % 'foilang_trend','foilang_seas','foilang_resid'
    flagtipY = (flagtipY_trend + flagtipY_resid)';
    foilang = (foilang_trend + foilang_resid)';
end

%% downsample to -- A (airfoil), F (flag)
time = time';
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

f1=figure;
subplot(3,1,1);
% Airfoil
plot(t(tidx),A(tidx), '.b','LineWidth',1.5,'MarkerSize',12); hold on;
% Plot raw ts for same time duration to compare
st_time = min(t(tidx)); end_time = max(t(tidx));
[st_time_A,st_idx_A] = min(abs(time-st_time));
[end_time_A,end_idx_A] = min(abs(time-end_time));
plot(time(st_idx_A:end_idx_A), foilang(st_idx_A:end_idx_A), ...
          '-k','LineWidth',1.5);
xlim([tminplot tmaxplot]); ylim([-50 50]);
ylabel('A');
% Flag
subplot(3,1,2); 
plot(t(tidx),F(tidx)*100, '.r','LineWidth',1.5,'MarkerSize',12); hold on; 
% Plot raw ts for same time duration to compare
plot(time(st_idx_A:end_idx_A), flagtipY(st_idx_A:end_idx_A)*100, ...
          '-k','LineWidth',1.5);
xlim([tminplot tmaxplot]); ylim([-3 3]); yticks([-3 0 3]);
ylabel('F_T');
% LDV
subplot(3,1,3); 
plot(t(tidx),u(tidx), '.m','LineWidth',1.5,'MarkerSize',12); hold on; 
% Plot raw ts of LDV for same time duration to compare
st_time = min(t(tidx)); end_time = max(t(tidx));
[st_time_LDV,st_idx_LDV] = min(abs(tmov-st_time));
[end_time_LDV,end_idx_LDV] = min(abs(tmov-end_time));
plot(tmov(st_idx_LDV:end_idx_LDV), umov(st_idx_LDV:end_idx_LDV), ...
          '-k','LineWidth',1.5);
xlim([tminplot tmaxplot]); ylim([0.1 0.7]); yticks([0.1 0.7]);
ylabel('u');
xlabel('time (s)');

f1.Position = [1 100 940 700];
saveas(gcf,[save_fig_dir filename '_processed_ts' '_ups' num2str(upsample)],'png');


%% Plot Raw time series
f0=figure;
subplot(3,1,1);
% Airfoil
st_time = min(t(tidx)); end_time = max(t(tidx));
[~,st_idx_A] = min(abs(time-st_time));
[~,end_idx_A] = min(abs(time-end_time));
plot(time(st_idx_A:end_idx_A), foilang_raw(st_idx_A:end_idx_A), ...
          '-k','LineWidth',1.5);
xlim([tminplot tmaxplot]); ylim([-50 50]);
ylabel('\theta (deg.)');
% Flag
subplot(3,1,2);
plot(time(st_idx_A:end_idx_A), flagtipY_raw(st_idx_A:end_idx_A)*100, ...
          '-k','LineWidth',1.5);
xlim([tminplot tmaxplot]); ylim([-3 3]); yticks([-3 0 3]);
ylabel('\Delta y_T (cm)');
% LDV
subplot(3,1,3); 
st_time = min(t(tidx)); end_time = max(t(tidx));
[~,st_idx_LDV] = min(abs(tLDV-st_time));
[~,end_idx_LDV] = min(abs(tLDV-end_time));
%plot(t(tidx),u(tidx), '-om','LineWidth',1.5,'MarkerSize',5); hold on; %title('LDV velocity (m/s)'); 
plot(tLDV(st_idx_LDV:end_idx_LDV), uLDV(st_idx_LDV:end_idx_LDV), ...
          '-k','LineWidth',1.5); 
xlim([tminplot tmaxplot]); ylim([0.1 0.7]); yticks([0.1 0.7]);
ylabel('U (m/s)');
xlabel('time (s)');

f0.Position = [1 100 940 700];
saveas(gcf,[save_fig_dir filename '_raw_ts' '_ups' num2str(upsample)],'png');


% Calculate the max and SD of deflection
disp("F_T: max = "+num2str(max(flagtipY_raw)) + ", SD = " + num2str(std(flagtipY_raw)))

%% Symbolize ts
piA = symbolize_data(A,m_embed);
piF = symbolize_data(F,m_embed);
piu = symbolize_data(u,m_embed);

% Plot symb time series

f8=figure;
subplot(3,1,1);
% Airfoil
plot(t(tidx),piA(tidx), '-b','LineWidth',1.5,'MarkerSize',12); hold on; %title('Airfoil angle (deg.)');
ylabel('\pi_A'); xlim([tminplot tmaxplot]); 
% Flag
subplot(3,1,2); 
plot(t(tidx),piF(tidx), '-r','LineWidth',1.5,'MarkerSize',12); hold on; %title('Flag deflection (cm)'); 
ylabel('\pi_{F_T}'); xlim([tminplot tmaxplot]); 
% LDV
subplot(3,1,3); 
plot(t(tidx),piu(tidx), '-m','LineWidth',1.5,'MarkerSize',12); hold on; %title('LDV velocity (m/s)'); 
ylabel('\pi_u'); xlim([tminplot tmaxplot]); 
xlabel('time (s)');

f8.Position = [1 100 940 700];
saveas(gcf,[save_fig_dir filename '_symb_ts' '_ups' num2str(upsample)],'png');


%% Transfer entropy

% Calc TE A->F and F->A symbolic
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
ylim([0 0.03]); yticks([0.0 0.01 0.02 0.03]);
f2.Position = [750 400 500 250];
saveas(gcf,[save_fig_dir filename '_delayTE' '_ups' num2str(upsample)],'png');

%% Pick the delay with max TE 
[maxTA_F,maxTA_F_idx] = max(TA_F);
del = delarr(maxTA_F_idx)
[maxnetTA_F,maxnetTA_F_idx] = max(netTA_F);
delB = delarr(maxnetTA_F_idx);
if (del~=delB)
    disp('Peaks of TE and net TE are at different delays');
end
TA_Fval = maxTA_F
TF_Aval = TF_A(maxTA_F_idx);

%% Shuffle source and test significance
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
xlabel('TE_{A \rightarrow {F_T}}'); ylabel('PDF'); %legend('','95 percentile','TE_{A->F}');
f4.Position = [350 400 400 250];
xlim([1e-5 1e-1]); xticks([1e-5 1e-4 1e-3 1e-2 1e-1]);
ylim([0 2500]); yticks([0 500 1000 1500 2000 2500]);
saveas(gcf,[save_fig_dir filename '_TE_AF' '_ups' num2str(upsample)],'png');

%% Cond TE_ A->F | u

TA_F_u = nan(Ndel,1);
for idel = 1:Ndel
    del2 = delarr(idel);
    TA_F_u(idel) = cond_transfer_entropy_delay(piF,piA,piu,del,del2,flag_norm); 
end

% Plot 
f6=figure;
plot(delarr(1:del),TA_Fval*ones(size(delarr(1:del),2),1),'--k','LineWidth',2); hold on;
plot(delarr(1:del),TA_F_u(1:del),'-r','LineWidth',2); 
xlabel('\delta_2'); %title('cond TE vs. delay of u wrt F');
ylabel('TE_{A \rightarrow {F_T}|u}');
ylim([0.02 0.024]); %yticks([0.040 0.041 0.042 0.043 0.044 0.045]);
f6.Position = [950 100 320 250];
saveas(gcf,[save_fig_dir filename '_delaycondTE_AFu' '_ups' num2str(upsample)],'png');

% Pick the delay with min cond TE 
[minTA_F_u,minTA_F_u_idx] = min(TA_F_u(1:del)); % Find del2 in (0,del), where del is max net TE A->F 
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
xlim([0.021 0.024]); xticks([0.021 0.022 0.023 0.024]);
ylim([0 4000]);
saveas(gcf,[save_fig_dir filename '_condTE_AFu' '_ups' num2str(upsample)],'png');