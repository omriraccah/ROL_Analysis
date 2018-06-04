function [ROL_peak,ROL_thresh]= Fit_ROL_TbT(path,elec,condition,deactive)

%% ROL Thresh

D = spm_eeg_load(path);

tt(1) = indsample(D,0/1000);
tt(2) = indsample(D,1500/1000);

idt = tt(1):tt(2);

active_means = [];
deactive_means = [];

thr_value = .5;


b_num = -.2;
b_value = indsample(D, b_num);
s_value = indsample(D, 0);
e_value = indsample(D, .8);

tr_toplot = setdiff(indtrial(D,condition),badtrials(D));

ROL_peaks = zeros(length(tr_toplot),length(elec));
ROL_thresh = zeros(length(tr_toplot),length(elec));


for e = 1:length(elec)

    %% ROL Peak

    for t = 1:length(tr_toplot)

        mdata = D(elec(e),idt,tr_toplot(t));
        baseline = mdata(b_value:s_value);
        b_mean = mean(baseline);
        b_std = std(baseline);
        signal = mdata(s_value:e_value);

        if deactive(e)

            % Peak
            pass = find(signal == min(signal(:)));
            ROL_peaks(t,e) = D.time(pass+s_value);
            
            % Threshold
            thr = b_mean-b_std*thr_value;
            pass = find(signal < thr);
            
            if isempty(pass)
                
                ROL_thr(t,e) = NaN;
            
            else
                
                ROL_thr(t,e) = D.time(pass(1)+s_value);
            
            end

        else

            % Peak
            pass = find(signal == max(signal(:)));
            ROL_peaks(t,e) = D.time(pass+s_value);
            
            % Threshold
            thr = b_mean+b_std*thr_value;
            pass = find(signal > thr);
            
            if isempty(pass)
                
                ROL_thr(t,e) = NaN;
            
            else
                
                ROL_thr(t,e) = D.time(pass(1)+s_value);
            
            end


        end


    end

end


ROL_peak = ROL_peaks;
ROL_thresh = ROL_thr;

