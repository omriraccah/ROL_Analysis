
chans = [107];

%chans = [90,73,89,72,88,71,87,70,86,69,85,68,84,67,83,66,82];

D = spm_eeg_load;

labs = 'autobio';
tt_b(1) = indsample(D,-.1);
tt_b(2) = indsample(D,0);

tt_s(1) = indsample(D,0);
tt_s(2) = indsample(D,1);

p_values = [];

for i = 1:length(chans)
 
s_values = [];
b_values = [];
    
    
    tr_toplot = setdiff(indtrial(D,labs),badtrials(D)); %take bad trials out
    %s_data = D(chans(i),tt_s(1):tt_s(2),tr_toplot);
    %b_data = D(chans(i),tt_b(1):tt_b(2),tr_toplot);

    
    % loop through trials
    for t = 1:length(tr_toplot)
    
        s_value = mean(D(chans(i),tt_s(1):tt_s(2),tr_toplot(t)));
        b_value = mean(D(chans(i),tt_b(1):tt_b(2),tr_toplot(t)));
        
        s_values = [s_values, s_value];
        b_values = [b_values, b_value];   
        
    end
    
    p_actual = permutation_paired(b_values, s_values,50000);
    
    p_values = [p_values,min(p_actual,1-p_actual)*2];

end



