function [onsets,peaks,slope_peaks] = Fit_ROL_TbT_Bins(path,elec,condition,deactive)

%% ROL Thresh

D = spm_eeg_load(path);

tt(1) = indsample(D,0/1000);
tt(2) = indsample(D,800/1000);

idt = tt(1):tt(2);

active_means = [];
deactive_means = [];


end_value = .5;
b_num = -.2;
b_value = indsample(D, b_num);
s_value = indsample(D, 0);
e_value = indsample(D, end_value);

tr_toplot = setdiff(indtrial(D,condition),badtrials(D));

ROL_peaks = zeros(length(tr_toplot),length(elec));
ROL_thresh = zeros(length(tr_toplot),length(elec));

params = struct; 
params.thr = 1;
params.thr_counter = 20;  
params.bin_timewin = 0.002;
times_nonover = 25;

onsets = zeros(length(tr_toplot),length(elec));
peaks = zeros(length(tr_toplot),length(elec));
slope_peaks = zeros(length(tr_toplot),length(elec));

%% Save Bin Data

for e = 1:length(elec)


    for t = 1:length(tr_toplot)
        
        mdata = D(elec(e),idt,tr_toplot(t));
        baseline = mdata(b_value:s_value);
        b_mean = mean(baseline);
        b_std = std(baseline);
        signal = mdata(s_value:e_value);
        

        % "params.bin_timewin" ms. wide bins
        numbins = (end_value / params.bin_timewin)-times_nonover;

        % Create cell structure to hold bin data
        bin_data = cell(floor(numbins),5); % bin_data will contain (data,averages, slopes,statistics), index of first value
        stim_onset_inx = indsample(D, 0);

        
        for b = 1:numbins
    
            if deactive(e) == 0


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
                if bin_data{b,2} > ((b_std*params.thr)+b_mean)

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
                if bin_data{b,2} < (b_mean-(b_std*params.thr))

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

        %% Get ROL

        ROL_Bin = (peak_bin_index-params.thr_counter)+1;

        if (peak_bin_index-params.thr_counter) < 0

            onset = NaN;
            
        else
            
            onset_index = bin_data{ROL_Bin+1,5};
            onset = D.time(onset_index);
            

        end
%         
        if onset < 0.05
%             
              onsets(t,e) = NaN;
         else

            onsets(t,e) = onset;
            
            
        end
        
        % get peaks
        if vdeactive == 0
            
         index = find(cell2mat(bin_data(:,2)) == max(cell2mat(bin_data(:,2))));
         onset_index = bin_data{index,5};
         onset = D.time(onset_index);   
            
         peaks(t,e) = onset;
         
        else
            
         index = find(cell2mat(bin_data(:,2)) == min(cell2mat(bin_data(:,2))));
         onset_index = bin_data{index,5};
         onset = D.time(onset_index);   
            
         peaks(t,e) = onset;
         
        end
        
        % get slope peaks
        if deactive == 0
            
         index = find(cell2mat(bin_data(:,3)) == max(cell2mat(bin_data(:,3))));
         onset_index = bin_data{index,5};
         onset = D.time(onset_index);   
            
         slope_peaks(t,e) = onset;
         
        else
            
         index = find(cell2mat(bin_data(:,3)) == min(cell2mat(bin_data(:,3))));
         onset_index = bin_data{index,5};
         onset = D.time(onset_index);   
            
         slope_peaks(t,e) = onset;
         
        end
        
        
    end
    
    
end



