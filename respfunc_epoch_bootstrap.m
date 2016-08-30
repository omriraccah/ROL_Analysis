function [Resp_data]= respfunc_epoch_bootstrap(D,elec,condition,params)

%% Script Info

% data: HFB timecourse for single electrode (all trials)
% onsets: stimulus onset times (s)
% durs: duration of each trial (s)
% maxdur: maximum trial length to consider (s)
% srate: sampling rate (Hz)
% preevent: time before onset to be used for baseline distribution (s)
% bs_frac: fraction of trials to be used in each boostrap repetition (0-1)
% bs_reps: # of repetitions for boostrapping
% params.smwin: gaussian smoothing window width (s)
% params.thr: threshold (# of SDs above mean of baseline data)
% params.thrwin: window size to use in sliding window (does average of 
%                   window surpass threshold?)
% params.mindur: average of thrwin must surpass threshold for at least 
%                   mindur (s) consecutive windows
% params.winlen: length of window for determining lines of best fit on data
%                   around threshold crossing (s)
% params.overlap: length of overlap between consencutive windows for line
%                   of best fit (s)

% Mo Dastjerdi - LBCN

% Updates for ver3, BLF - LBCN 2014.
% Cleaned up code, added comments
% Modified moving window, larger with more overlap

% modified by Omri Raccah- 2016

%% set default values for inputs (if needed)
if isempty(D)
    D = spm_eeg_load;
end



% if exist(file_name) == 0
%     
%     D = spm_eeg_load;
%     
% elseif exist(file_name) == 1
%     
%     D = spm_eeg_load(file_name)
%     
% end


% if nargin<4
%     srate = 1000;
% end
srate = D.fsample;

% if nargin<5
%     bs_frac = 0.5;
% end

% if nargin<4 
%     bs_reps = 1000;
% end

%% set default values for params if undefined

if nargin<4
    
    params = struct;
    
    if ~isfield(params,'smwin')
        params.smwin = 0.005;   
    end
    if ~isfield(params,'thr')
        params.thr = 5;         
    end
    if ~isfield(params,'thrwin')
        params.thrwin = 0.05;    
    end
    if ~isfield(params,'mindur')
        params.mindur = 0.05;
    end
    if ~isfield(params,'winlen')
        params.winlen = 0.07;
    end
    if ~isfield(params,'overlap')
        params.overlap = 0.065;
    end
 
end

% make figure directory
mkdir('figs')

% Loop all code through electrodes
% for e = 1:length(elec)

% isolate codition data
data= squeeze(D(elec,:,condition));
data=data';
time=D.time; 


%%
% onset_pt = floor(-preevent*srate);
maxdur= time(end);
maxdur_pt = floor(maxdur*srate);
ntrials = size(data,1);
bastw = 1:indsample(D,0);

thrwin_pt = floor(params.thrwin*srate);
mindur_pt = floor(params.mindur*srate);

%% smooth data with gaussian
%
% winSize = floor(srate*params.smwin);
% gusWin= gausswin(winSize)/sum(gausswin(winSize));
% sm_data = convn(data,gusWin','same');


%% determine threshold based on null/baseline (ITI) distribution

% baseline_bs = []; %bootstrapped baseline
%
% ntrials_bs = ceil(ntrials*bs_frac);
%
% bs_trials = zeros(bs_reps,ntrials_bs);
%
% for ri = 1:bs_reps
%     bs_trials(ri,:)=randperm(ntrials,ntrials_bs);
% end

baseline = data(:,bastw);



% for ri = 1:bs_reps
%     temp_trials = bs_trials(ri,:);
%     baseline_bs = [baseline_bs nanmean(baseline(temp_trials,:))];
% end

mean_bl = nanmean(baseline(:));
std_bl = nanstd(nanstd(baseline,0,1));

% mean_bs = nanmean(baseline_bs(:));
% std_bs = nanstd(baseline_bs(:));

thr_val = mean_bl + params.thr*std_bl;

%baseline correct
trial_data_bc = (data-mean_bl)/std_bl;

%% reorganize data into trials (include ITI)
% trial_data = NaN*ones(ntrials,maxdur_pt+length(bastw)+1);
% trial_data_bc = NaN*ones(ntrials,maxdur_pt+length(bastw)+1);
% for ti = 1:ntrials
%     if (durs_pt(ti)>maxdur_pt)
%         trial_data(ti,:)=sm_data(onset_pt(ti)-length(bastw):onset_pt(ti)+maxdur_pt);
%         trial_data_bc(ti,:)=sm_data_bc(onset_pt(ti)-length(bastw):onset_pt(ti)+maxdur_pt);
%     else
%         trial_data(ti,1:length(bastw)+durs_pt(ti)+1) = sm_data(onset_pt(ti)-length(bastw):onset_pt(ti)+durs_pt(ti));
%         trial_data_bc(ti,1:length(bastw)+durs_pt(ti)+1) = sm_data_bc(onset_pt(ti)-length(bastw):onset_pt(ti)+durs_pt(ti));
%     end
% end

%% jackknifing
bef_time = 0.2; %time before/after threshold crossing pt to use for line-fitting
aft_time = 0.1;

% parameters for line fitting
% winlen = 0.1; % in seconds
% overlap = 0.09;

% winlen = 0.07; % in seconds
% overlap = 0.06;

bef_pt = floor(bef_time*srate);
aft_pt = floor(aft_time*srate);

winlen_pt = floor(params.winlen*srate);
overlap_pt = floor(params.overlap*srate);

% bs_onsets = NaN*ones(1,bs_reps); %bootstrapped (jackknifed) onsets
% bs_onsets = zeros(1,bs_reps);

% trial_data_bs = NaN*ones(bs_reps,length(time));

% if (size(data,1)>1)
%     for ti = 1:bs_reps
%         temp_trials = bs_trials(ti,:);
temp_mn = nanmean(data);
temp_mn = (temp_mn-mean_bl)/std_bl;
%          trial_data_bs(ti,:)=temp_mn;
numwins = 0; % number of consecutive windows surpassing thresh
ind = 1; %
bs_onsets = NaN;
while ((numwins < mindur_pt) && ind<(maxdur_pt+length(bastw)-thrwin_pt+1))
    temp_inds = ind:ind+thrwin_pt-1;
    if (nanmean(temp_mn(temp_inds))>params.thr)
        %             if (nanmean(temp_mn(temp_inds))>thr_val)
        numwins = numwins+1;
    else
        numwins = 0;
    end
    ind = ind+1;
end
if (numwins >= mindur_pt)
    ind = ind-mindur_pt; % first point that threshold was crossed
    if (ind>bef_pt)
        if (ind<maxdur_pt+length(bastw)-aft_pt)
            win_pt = ind-bef_pt:ind+aft_pt;
        else
            win_pt = ind-bef_pt:maxdur_pt+length(bastw)+1;
        end
    else
        win_pt = 1:ind+aft_pt;
    end
    
    % split into overlapping windows, find line of best fit for each segment
    sig_tmp = temp_mn(win_pt);
    w_time = time(win_pt);
    sig_tmp= buffer(sig_tmp,winlen_pt,overlap_pt,'nodelay');
    t_tmp= buffer(w_time,winlen_pt,overlap_pt,'nodelay');
    sig_tmp(:,end)=[];
    t_tmp(:,end)=[];
    nwins = size(sig_tmp,2);
    
    slopes=NaN*ones(1,nwins);
    mse=NaN*ones(1,nwins);
    for ii=1:nwins
        Ps= polyfit(t_tmp(:,ii),sig_tmp(:,ii),1);
        y2= polyval(Ps,t_tmp(:,ii));
        slopes(ii)= Ps(1);
        mse(ii)= sum((sig_tmp(:,ii)-y2).^2);
    end
    
    [s_tmp iA]= sort(slopes,'descend'); %slopes
    [e_tmp iB]= sort(mse,'ascend'); %errors
    %get smallest error, of top 5 slopes
    i_tmp= find( mse(iA(1:5))== min(mse(iA(1:5))) );
    %         i_tmp = find(slopes(iB(1:5))==max(slopes(iB(1:5))));
    if ~isempty(i_tmp)
        bs_onsets= t_tmp(1,iA(i_tmp));
    end
end
%     end
% end


%% Output data
Resp_data.all = bs_onsets;
% Resp_data.mn = nanmean(bs_onsets);
% Resp_data.sd = nanstd(bs_onsets);
Resp_data.trace_bc = trial_data_bc;
% Resp_data.trace_bs = trial_data_bs;
% Resp_data.meantrace_bs = (nanmean(data)-mean_bs)/std_bs;
Resp_data.time = time;

%% plot trial averaged signal
% close all
% plot(time,nanmean(data),'k-','LineWidth',3)
% hold on
% plot([Resp_data.mn,Resp_data.mn],ylim,'r-','LineWidth',2)
%
% %Response onset time for all trials in secs
% Resp_data.all_onsets_sec = rsp_onset;
% %HFB value for onset all trials
% Resp_data.all_onset_val_HFB = rsp_val;
% %Max HFB value all trials
% Resp_data.all_max_resp = max_rsp';
% %Time of max HFB value all trials
% Resp_data.all_max_resp_time_sec = max_rsp_time';
% %Duration of HFB above threshold (total) in sec
% Resp_data.all_tot_resp_dur_sec = tot_rsp_dur';
% %Duration of HFB above threshold (total) in percent of trial duration
% Resp_data.tot_resp_dur_prct = tot_rsp_dur_prct';
% %Cummulative response curves of HFB all trials
% Resp_data.all_cum_resp = cum_rsp';
% %Duration of each trial in sample points
% Resp_data.all_trial_length_pts = trial_durs;
% %sample rate of data
% Resp_data.srate = srate;

%% Remove Outlier
% 
% % 
% all_resp = Resp_data.all;
% 
% 
% noise_std = 2;
% 
% pre_outlier= zscore(all_resp)
% 
% a= find (pre_outlier > noise_std | pre_outlier < -(noise_std))
% 
% % a = find (pre_outlier) >= pre_outlier + pre_outlier | pre_outlier - noise_std >= pre_outlier
% all_resp(a) = nan;


% 
% rsp_outliers = find (all_resp) >= (nanmean(abs(all_resp))) + noise_std*nanstd (abs(all_resp))|(nanmean(abs(all_resp)))...
%     - noise_std*nanstd (abs(all_resp)) >= all_resp;
% 




%% remove outlier; change the noise_std = 2; add it to the function


% end
