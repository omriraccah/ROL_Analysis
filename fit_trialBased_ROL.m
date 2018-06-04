function [onsets,peaks,slope_peaks] = fit_trialBased_ROL(path,elec,condition,deactive,norm)
% edited by Omri on Mar. 23, 2018

%% ROL thr_valueesh

% Load data
D = spm_eeg_load(path);

% Get "good" trials
tr_toplot = setdiff(indtrial(D,condition),badtrials(D));

% Set ROL parameters
params = struct; 
params.thr_value = 1; % How many SD over baseline period
params.thr_value_counter = 25; % How many bins must surpass threshold
params.bin_timewin = 0.002; % How much non-overlap between bins
times_nonover = 15; % 30ms bins with 28ms overlap

% Set range parameters
end_value = .5; % 500 ms
start_value = 0; % stim_onset
bas_value = 0; % end of baseline
stim_onset_inx = indsample(D, start_value);
stim_bas_inx = indsample(D, bas_value);
b_start = indsample(D, -0.2);
end_time_inx = indsample(D, end_value);

% Set data structures
onsets = NaN*zeros(length(tr_toplot),length(elec));
peaks = NaN*zeros(length(tr_toplot),length(elec));
slope_peaks = NaN*zeros(length(tr_toplot),length(elec));

%% Calculate Bin Data

% Loop through electrodes
for e = 1:length(elec)
    
    % Load data for specific electrode and trial
    alldata = squeeze(D(elec(e),:,tr_toplot));
    
    % Normalize trial
    win_max = stim_onset_inx:end_time_inx;
    m = max(alldata(win_max,:),[],1);
    if norm
        alldata = alldata./repmat(m,size(alldata,1),1);
    end
    
    % Calculate baseline mean and standard deviation
    baseline_signal = mean(alldata(b_start:stim_bas_inx,:),2);
    mean_bl = mean(baseline_signal);
    std_bl = std(baseline_signal);

    % Loop through trials
    for t = 1:length(tr_toplot)    
        
        % Load data for specific trial
        mdata = alldata(:,t)';
                
        % Get number of bins
        numbins = (((end_value-start_value) / params.bin_timewin)-times_nonover)+1;

        % Create cell structure to hold bin data
        bin_data = cell(floor(numbins),6); % bin_data will contain (data; averages; slopes; whether average surpasses threshold; index of first bin value; index of last bin value) 

        % Loop through bins and create data structure
        for b = 1:numbins
    
            % IF ACTIVE SITE (as appose to deactive)
            if deactive(e) == 0

                % Set threshold value
                thr = mean_bl + params.thr_value*std_bl;

                if b == 1
            
                    % Set data segment of interest and store signal data
                    s_value = indsample(D, start_value);
                    e_value = indsample(D, start_value+(times_nonover*params.bin_timewin));
                    bin_data(b,1) = {mdata((s_value):((e_value)))};

                else

                    % Set data segment of interest and store signal data
                    s_value = indsample(D, start_value+ (params.bin_timewin*(b-1)));
                    e_value = indsample(D, start_value+(params.bin_timewin*(b-1))+(times_nonover*params.bin_timewin));
                    bin_data(b,1) = {mdata(((s_value):((e_value))))};
                    
                   % one can run the line below to ensure that the bin sizes are generating correctly through the loop 
                   % bin_range = [start_value+ (params.bin_timewin*(b-1)),start_value+(params.bin_timewin*(b-1))+(times_nonover*params.bin_timewin)]
                   
                end

                % Create cell structure and store bin average
                bin_data(b,2) = {mean(bin_data{b,1})};

                % Get slope
                coefficients = polyfit(1:length(bin_data{b,1}), bin_data{b,1}, 1);
                slope = coefficients(1);
                bin_data(b,3) = {slope};

                % Check if bin mean exceeds thershold (1 = yes, 0 = no)
                if bin_data{b,2} > thr

                    bin_data(b,4) = {1};

                else

                    bin_data(b,4) = {0};

                end

                % Get index of first value
                bin_data(b,5) = {s_value};
                
                % Get index of last value
                bin_data(b,6) = {e_value};
                
                
            % IF DEACTIVE SITE  
            else
        
                % Set threshold value
                thr = mean_bl - params.thr_value*std_bl;
                
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

                % Check if bin mean exceeds thr_valueeshold (1 = yes, 0 = no)
                if bin_data{b,2} < thr

                    bin_data(b,4) = {1};

                else

                    bin_data(b,4) = {0};

                end

                % Get index of first value
                bin_data(b,5) = {s_value};
                
                % Get index of last value
                bin_data(b,6) = {e_value};
          
                
            end
            
        end

        %% Estimate ROL

        % counter to track slope and threshold 
        counter = 0;
        peak_bin_index = 0;

        % loop through bin data to find  ten consequtive bins meeting threshold
        % criterion
        for i = 1:length(bin_data) % would be better with 'while'

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

           
            end

        end

        
        %% Initialize ROL

        ROL_Bin = (peak_bin_index-params.thr_value_counter)+1;

        if (peak_bin_index-params.thr_value_counter) >= 0
            
            onset_index = bin_data{ROL_Bin,5};
            time = D.time;
            onset = time(onset_index);
            
            % Estimate signal's peak from bins
            % win_peak = ROL_Bin:size(bin_data,1);
            if ~deactive
                [~, Peak_Bin]=max([bin_data{:,2}]);
                peak_index = bin_data{Peak_Bin,5};
                peak = time(peak_index);
                
            else
                [~, Peak_Bin]=min([bin_data{:,2}]);
                peak_index = bin_data{Peak_Bin,5};
                peak = time(peak_index);
            end
            peaks(t,e) = peak;                                   
            onsets(t,e) = onset;
            
            % show trial-by-trial onset
            
            % Un-comment contents to see trial by trial estimates
            if ~isnan(onset) 
               %figure; plot(D.time(stim_onset_inx:end_time_inx), mdata(stim_onset_inx:end_time_inx));
               %hold on; plot(onset, .3,'r*') 
               %hold on; plot(peak, .3,'g*') 
                
            end
            
            % Estimate slope in window of positivity to reach onset           
            %twin = bin_data{ROL_Bin+tpeak-1,6};
            %slope_peaks(t,e) = (mdata(twin)-mdata(onset_index))/(length(onset_index:twin)*(1/D.fsample));
            
            
        end

        
    end
    
    
end
