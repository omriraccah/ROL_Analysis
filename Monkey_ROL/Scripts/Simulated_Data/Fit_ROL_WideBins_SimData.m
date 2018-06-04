function [ROL_Data]= Fit_ROL_WideBins_SimData(D,elec,condition,params,deactive,mean_signal)
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
b_num = 0;
b_value = indsample(D, b_num);
s_value = indsample(D, .1);

    params = struct;
    
    if ~isfield(params,'bastw')
        params.bastw = b_value:s_value;         
    end
    if ~isfield(params,'thr')
        params.thr = 2;         
    end
    if ~isfield(params,'thr_counter')
        params.thr_counter = 30;         
    end
    if ~isfield(params,'thrwin')
        params.thrwin = 0.05;    
    end
    if ~isfield(params,'mindur')
        params.mindur = 0.05;
    end
    if ~isfield(params,'bin_timewin')
        params.bin_timewin = 0.002;
    end
 
end_value = .5;
    
% Set data segment of interest
e_value = indsample(D, end_value);

% Get signal data
data= squeeze(D(elec,:,condition));
% Get time
time=D.time(s_value:e_value);

%% Initialize variables

data = data';

% Get baseline values, average, and standard deviation
baseline = data(:,params.bastw);
mean_bl = nanmean(baseline(:));
mean_across_trials = nanmean(baseline,1);
std_bl = nanstd(mean_across_trials);

% Define threshold value
thr_val = mean_bl + params.thr*std_bl;

% Average data across trials
mdata = nanmean(data);
mdata = mdata(s_value:e_value);

%baseline correct
mdata = (mdata-mean_bl)/std_bl;


%% Create Data structure for individual bins

times_nonover = 30;

% "params.bin_timewin" ms. wide bins
numbins = (end_value / params.bin_timewin)-times_nonover;

% Create cell structure to hold bin data
bin_data = cell(floor(numbins),5); % bin_data will contain (data,averages, slopes,statistics), index of first value
stim_onset_inx = indsample(D, .1);

% Save bin data

for b = 1:numbins
    
    if deactive == 0
        
        if b == 1
            
            % Set data segment of interest and grab values
            s_value = indsample(D, .1);
            e_value = indsample(D, params.bin_timewin+(times_nonover*params.bin_timewin));
            bin_data(b,1) = {mdata(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
        else
            
            % Set data segment of interest and grab values
            s_value = indsample(D, params.bin_timewin*(b-1))+1;
            e_value = indsample(D, params.bin_timewin*(b)+(times_nonover*params.bin_timewin));
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
            s_value = indsample(D, .1);
            e_value = indsample(D, params.bin_timewin+(times_nonover*params.bin_timewin));
            bin_data(b,1) = {mdata(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
        else
            
            % Set data segment of interest and grab values
            s_value = indsample(D, params.bin_timewin*(b-1))+1;
            e_value = indsample(D, params.bin_timewin*(b)+(times_nonover*params.bin_timewin));
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
        if bin_data{i,4} == 1 && bin_data{i,3} > 0

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
        if bin_data{i,4} == 1 && bin_data{i,3} < 0

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


%% NEW: GET BASELINE SLOPE BINS

% "params.bin_timewin" ms. wide bins
numbins_bl = (abs(b_num) / params.bin_timewin)-times_nonover;

% Create cell structure to hold bl_slope data
bin_data_bl = cell(floor(numbins_bl),2); % bin_data will contain (data,averages, slopes,statistics), index of first value
stim_onset_inx = indsample(D, b_num);

% Save bin data
for b = 1:numbins_bl
    
    if b == 1
            
            % Set data segment of interest and grab values
            s_value = indsample(D, b_num);
            e_value = indsample(D, b_num+params.bin_timewin+(times_nonover*params.bin_timewin));
            bin_data_bl(b,1) = {mean_across_trials(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
    else
            
            % Set data segment of interest and grab values
            s_value = indsample(D, b_num+params.bin_timewin*(b-1))+1;
            e_value = indsample(D, b_num+params.bin_timewin*(b)+(times_nonover*params.bin_timewin));
            bin_data_bl(b,1) = {mean_across_trials(((s_value-stim_onset_inx)+1):((e_value-stim_onset_inx)+1))};
            
    end
    
    % Get slope
    coefficients = polyfit(1:length(bin_data_bl{b,1}), bin_data_bl{b,1}, 1);
    slope = coefficients(1);
    bin_data_bl(b,2) = {slope};
   
end

%% Get bin index for ROLUnfor
ROL_Bin = (peak_bin_index-params.thr_counter)+1;
if (peak_bin_index-params.thr_counter) < 0
    
    onset = NaN;
    ROL_Data.all = onset;
    
    if mean_signal == 1
        
         ROL_Data.bin_data = {bin_data}; 
        
    end
    
    
else
 
    % FIND SLOPE THRESHOLD
    if deactive == 0

        % this should be a measure of time
        %M = max(mdata)/bin_timewin;
        %[sortedValues,sortIndex] = sort(cell2mat(bin_data(:,3)),'descend');  %# Sort the values in
        %SlopeValues = sortedValues(1:10);

        SlopeValues = cell2mat(bin_data(peak_bin_index-(params.thr_counter-1):peak_bin_index,3));
        slope_mean_c = mean(SlopeValues);
        slope_mean_b = mean(cell2mat(bin_data_bl(:,2)));
        slope_thr = ((slope_mean_c + slope_mean_b) / 2);
% %         
%         slope_std = std(cell2mat(bin_data_bl(:,2)));
%         slope_mean = mean(cell2mat(bin_data_bl(:,2)));
%         slope_thr = slope_mean + (slope_std*15);
%    
    else
        
        SlopeValues = cell2mat(bin_data(peak_bin_index-(params.thr_counter-1):peak_bin_index,3));
        slope_mean_c = mean(SlopeValues);
        slope_mean_b = mean(cell2mat(bin_data_bl(:,2)));
        slope_thr = ((slope_mean_c + slope_mean_b) / 2);
% 
%         slope_std = std(cell2mat(bin_data_bl(:,2)));
%         slope_mean = mean(cell2mat(bin_data_bl(:,2)));
%         slope_thr = slope_mean - (slope_std*15);
    
    end

    
    % onset_index is the first index in a bin -> used to grab the ROL value
    % from D.Time
 
    

    %run backword through bins to look at slope value
    for i = 1:ROL_Bin-1

          if deactive == 0

                if bin_data{ROL_Bin-i,3} < slope_thr


                    ROL_Bin = (ROL_Bin-i)+1;
                    break


                end

          end

          if deactive == 1

                if bin_data{ROL_Bin-i,3} > slope_thr  


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

