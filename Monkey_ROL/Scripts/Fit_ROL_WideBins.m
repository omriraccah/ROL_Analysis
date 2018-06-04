function [ROL_Data]= Fit_ROL_WideBins(D,elec,condition,params,deactive,mean_signal)
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



if deactive == 0
                
     params.thr = (mean_across_trials + params.thr*std_bl);
    
else
             
    params.thr = (mean_across_trials - params.thr*std_bl);
    
end
    

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
        
        % Get slope
        coefficients = polyfit(1:length(bin_data{b,1}), bin_data{b,1}, 1);
        slope = coefficients(1);
        bin_data(b,3) = {slope};

        % Check if bin mean exceeds threshold (1 = yes, 0 = no)
        if bin_data{b,2} > params.thr %&& bin_data{b,3} > 0
            
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
        
        % Get slope
        coefficients = polyfit(1:length(bin_data{b,1}), bin_data{b,1}, 1);
        slope = coefficients(1);
        bin_data(b,3) = {slope};
        
        % Check if bin mean exceeds threshold (1 = yes, 0 = no)
        if bin_data{b,2} < params.thr %&& bin_data{b,3} < 0 %mean(mean_across_trials)-(std_bl*params.thr)
            
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
    if ~deactive
        
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


%% NEW: GET BASELINE SLOPE BINS

e_value = indsample(D, 0 );
s_value = indsample(D, b_num);
coefficients = polyfit(1:length(mdata(s_value:e_value)), mdata(s_value:e_value), 1);
baseline_slope = coefficients(1);

%% Get bin index for ROLUnfor
ROL_Bin = (peak_bin_index-params.thr_counter)+1;
if (peak_bin_index-params.thr_counter) < 0
    
    onset = NaN;
    ROL_Data.all = onset;
    
    if mean_signal == 1
        
         ROL_Data.bin_data = {bin_data}; 
        
    end
    
    
else
 
    % Define slope threshold
    if ~deactive

        SlopeValues = cell2mat(bin_data(ROL_Bin:peak_bin_index,3));
        [~, max_bin] = max(SlopeValues);
        slope_mean_c = max(SlopeValues);
        slope_mean_b = baseline_slope; 
        slope_thr = ((slope_mean_c + slope_mean_b) / 2);
        % Reset ROL_Bin to steepest slope
        ROL_Bin = ROL_Bin+(max_bin-1);

    else
        
        SlopeValues = cell2mat(bin_data(ROL_Bin:peak_bin_index,3));
        [~, min_bin] = min(SlopeValues);
        slope_mean_c = min(SlopeValues);
        slope_mean_b = baseline_slope; 
        slope_thr = ((slope_mean_c + slope_mean_b) / 2);
        % Reset ROL_Bin to steepest slope
        ROL_Bin = ROL_Bin+(min_bin-1);
    
    end

    
    %run backwards through bins to look at slope value
    for i = 1:(ROL_Bin-1)

         % if activated elec
          if ~deactive

                if bin_data{ROL_Bin-i,3} < slope_thr


                    ROL_Bin = (ROL_Bin-i);
                    break


                end

          end

          if deactive

                if bin_data{ROL_Bin-i,3} > slope_thr 


                    ROL_Bin = (ROL_Bin-i);
                    break


                end

          end


    end


    % Get onset time
    onset_index = bin_data{ROL_Bin,5};
    onset = D.time(onset_index);
    ROL_Data.all = onset;


end

