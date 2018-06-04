function [ROL_Data]= Fit_ROL(D,elec,condition,params,deactive,mean_signal)
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
s_value = indsample(D, 0);

s_value = indsample(D, 0);

    params = struct;
    
    if ~isfield(params,'bastw')
        params.bastw = 1:s_value;         
    end
    if ~isfield(params,'thr')
        params.thr = 1;         
    end
    if ~isfield(params,'thr_counter')
        params.thr_counter =4;         
    end
    if ~isfield(params,'thrwin')
        params.thrwin = 0.05;    
    end
    if ~isfield(params,'mindur')
        params.mindur = 0.05;
    end
    if ~isfield(params,'bin_timewin')
        params.bin_timewin = 0.015;
    end
 
end_value = .8;
    
% Set data segment of interest
s_value = indsample(D, 0);
e_value = indsample(D, end_value);

% Get signal data
data= squeeze(D(elec,:,condition));
% Get time
time=D.time(s_value:e_value);

%% Normaliza data

% %% Normalize trial-by-trial
% % 
% if deactive == 0
%     
%     for n = 1:length(data(:,1))
% 
%         M = max(data(n,:));
%         data(n,:) = data(n,:)/M;
% 
%     end
%     
% end



%     
% elseif deactive == 1
%    
%       for n = 1:length(data(:,1))
% 
%         M = min(data(n,:));
%         data(n,:) = data(n,:)/abs(M);
% 
%       end
%     
% end
%% Initialize variables

data = data';

% Get baseline values, average, and standard deviation
baseline = data(:,params.bastw);
mean_bl = nanmean(baseline(:));
std_bl = nanstd(baseline(:));

% Define threshold value
thr_val = mean_bl + params.thr*std_bl;

% Average data across trials
mdata = nanmean(data);
mdata = mdata(s_value:e_value);

%baseline correct
mdata = (mdata-mean_bl)/std_bl;


%% Create Data structure for individual bins

% "params.bin_timewin" ms. wide bins
numbins = (end_value / params.bin_timewin);

% Create cell structure to hold bin data
bin_data = cell(floor(numbins),5); % bin_data will contain (data,averages, slopes,statistics), index of first value
stim_onset_inx = indsample(D, 0);

% Save bin data

for b = 1:numbins
    
    if deactive == 0
        
        if b == 1
            
            % Set data segment of interest and grab values
            s_value = indsample(D, 0);
            e_value = indsample(D, params.bin_timewin);
            bin_data(b,1) = {mdata(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
        else
            
            % Set data segment of interest and grab values
            s_value = indsample(D, params.bin_timewin*(b-1))+1;
            e_value = indsample(D, params.bin_timewin*(b));
            bin_data(b,1) = {mdata(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
        end
        
        % Create cell structure to hold bin averages
        bin_data(b,2) = {mean(bin_data{b,1})};
        
        % Get slope
        coefficients = polyfit(1:length(bin_data{b,1}), bin_data{b,1}, 1);
        slope = coefficients(1);
        bin_data(b,3) = {slope};

        % Check if bin mean exceeds threshold (1 = yes, 0 = no)
        if bin_data{b,2} > ((std_bl*params.thr)+mean_bl)
            
            bin_data(b,4) = {1};
            
        else
            
            bin_data(b,4) = {0};
            
        end
        
        % Get index of first value
        bin_data(b,5) = {s_value};
            
        
    else
        
        if b == 1
            % Set data segment of interest and grab values
            s_value = indsample(D, 0);
            e_value = indsample(D, params.bin_timewin);
            bin_data(b,1) = {mdata(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
        else
            
            % Set data segment of interest and grab values
            s_value = indsample(D, params.bin_timewin*(b-1))+1;
            e_value = indsample(D, params.bin_timewin*(b));
            bin_data(b,1) = {mdata(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
        end
        
        % Create cell structure to hold bin averages
        bin_data(b,2) = {mean(bin_data{b,1})};
        
        % Get slope
        coefficients = polyfit(1:length(bin_data{b,1}), bin_data{b,1}, 1);
        slope = coefficients(1);
        bin_data(b,3) = {slope};
        
        % Check if bin mean exceeds threshold (1 = yes, 0 = no)
        if bin_data{b,2} < (mean_bl-(std_bl*params.thr))
            
            bin_data(b,4) = {1};
            
        else
            
            bin_data(b,4) = {0};
            
        end
        
        % Get index of first value
        bin_data(b,5) = {s_value};
        
    end
    
    
end

%% Estimate ROL

% counter to track slope and threshold 
counter = 0;
peak_bin_index = 0;

% loop through bin data to find  ten consequtive bins meeting threshold
% criterion
for i = 1:length(bin_data)
    
    % if active electrode
    if deactive == 0
        
        % Check if bin is above threshold, with a positive slope
        if bin_data{i,4} == 1 

            counter = counter + 1;

            % Check if counter exceeds window threshold
            if counter ==  params.thr_counter

                % save the value of tenth bin
                peak_bin_index = i;

                break

            end


        else

            counter = 0;

        end
 
    % if deactive electrode
    else
        
        % Check if bin is above threshold, with a negative slope
        if bin_data{i,4} == 1 

            counter = counter + 1;

            % Check if counter exceeds window threshold
            if counter ==  params.thr_counter

                % save the value of tenth bin
                peak_bin_index = i;

                break

            end

        else

            counter = 0;

        end
    end
    
end

% Get bin index for ROL
ROL_Bin = (peak_bin_index-params.thr_counter)+1;
if (peak_bin_index-params.thr_counter) < 0
    
    onset = NaN;
    ROL_Data.all = onset;
    
    if mean_signal == 1
        
         ROL_Data.bin_data = {bin_data}; 
        
    end
    
    
else
 
    % onset_index is the firt index in a bin -> used to grab the ROL value
    % from D.Time
    
    % run backword through bins to look at slope value
    for i = 1:ROL_Bin-1
        
      if deactive == 0
          
            if bin_data{ROL_Bin-i,2} <= (mean_bl+(std_bl*.3))


                ROL_Bin = (ROL_Bin-i)+1;
                break


            end
            
      end
      
      if deactive == 1
          
            if bin_data{ROL_Bin-i,2} >= (mean_bl+(std_bl*.3))


                ROL_Bin = (ROL_Bin-i)+1;
                break


            end
            
      end
        
        
    end
    
onset_index = bin_data{ROL_Bin,5};
    
onset = D.time(onset_index);
ROL_Data.all = onset;
    
    if mean_signal == 1
        
         ROL_Data.bin_data = bin_data; 
         ROL_Data.ROL_Bin_Index = {ROL_Bin}; 
         ROL_Data.pos_thr = (std_bl*params.thr)+mean_bl; 
         ROL_Data.neg_thr = mean_bl-(std_bl*params.thr);
        
    end

end

