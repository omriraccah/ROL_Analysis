function [T2P]= Fit_T2P_WideBins(D,elec,condition,params,deactive,mean_signal)
%% Parameter Information

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

% Omri Raccah - LBCN

% modified by Omri Raccah- 2016

%% set default values for inputs (if needed)
if isempty(D)
    D = spm_eeg_load;
end
srate = D.fsample;

%% set default values for params if undefined

% find 0 index
b_num = -.2;
b_value = indsample(D, b_num);
s_value = indsample(D, 0);

    params = struct;
    
    if ~isfield(params,'bastw')
        params.bastw = b_value:s_value;         
    end
    if ~isfield(params,'thr')
        params.thr = 1;         
    end
    if ~isfield(params,'thr_counter')
        params.thr_counter = 50;         
    end
    if ~isfield(params,'bin_timewin')
        params.bin_timewin = 0.002;
    end
 
end_value = 1;
    
% Set data segment of interest
e_value = indsample(D, end_value);

% Get signal data
data= squeeze(D(elec,:,condition));  

%% Initialize variables

data = data';

% Get baseline values, average, and standard deviation
baseline = mean(data(:,params.bastw));
mean_across_trials = nanmean(baseline(:));
std_bl = nanstd(baseline(:)); 

mdata = nanmean(data);



%% Create Data structure for individual bins

times_nonover = 30;

% "params.bin_timewin" ms. wide bins
numbins = (end_value / params.bin_timewin)-times_nonover;

% Create cell structure to hold bin data
bin_data = cell(floor(numbins),5); % bin_data will contain (data,averages, slopes,statistics), index of first value
stim_onset_inx = indsample(D, 0);

% Save bin data

for b = 1:numbins
    
    if deactive == 0
        
        if b == 1
            
                    % Set data segment of interest and grab values
                    s_value = indsample(D, 0);
                    e_value = indsample(D, (times_nonover*params.bin_timewin));
                    bin_data(b,1) = {mdata((s_value):((e_value)))};

       else

                    % Set data segment of interest and grab values
                    s_value = indsample(D, params.bin_timewin*(b-1));
                    e_value = indsample(D, (params.bin_timewin*(b-1))+(times_nonover*params.bin_timewin));
                    bin_data(b,1) = {mdata(((s_value):((e_value))))};

                
        end
        
        % Create cell structure to hold bin averages
        bin_data(b,2) = {mean(bin_data{b,1})};
        
        
        % Get index of first value
        bin_data(b,3) = {s_value};
            
        
    else
        
        if b == 1
            
                    % Set data segment of interest and grab values
                    s_value = indsample(D, 0);
                    e_value = indsample(D, (times_nonover*params.bin_timewin));
                    bin_data(b,1) = {mdata((s_value):((e_value)))};

       else

                    % Set data segment of interest and grab values
                    s_value = indsample(D, params.bin_timewin*(b-1));
                    e_value = indsample(D, (params.bin_timewin*(b-1))+(times_nonover*params.bin_timewin));
                    bin_data(b,1) = {mdata(((s_value):((e_value))))};

                
        end
        
        
        % Create cell structure to hold bin averages
        bin_data(b,2) = {mean(bin_data{b,1})};
        
        
        % Get index of first value
        bin_data(b,3) = {s_value};
        
    end
    
    
end

%% Estimate ROL

if deactive == 0
    
    [~,m_inx] = max(cell2mat(bin_data(:,2)));
    s_inx = bin_data{m_inx,3};
    T2P = D.time(s_inx);
    
    
else
    
    
    [~,m_inx] = min(cell2mat(bin_data(:,2)));
    s_inx = bin_data{m_inx,3};
    T2P = D.time(s_inx);
    
    
end

