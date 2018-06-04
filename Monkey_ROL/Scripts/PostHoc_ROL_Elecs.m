elecs = indchan;

elec_size = size(onsets);


take_out = cell(1,1);
take_out_ind = 1;

keep_in = cell(1,1);
keep_in_ind = 1;

ROL_Values = [];


% loop through electrodes
for e = 1:elec_size(2)
    
    nan_counter = 0;
    
    for b = 1:elec_size(1)
    
        if isnan(onsets(b,e))
        
            nan_counter = nan_counter+1;
            
        end
        
    end
    
    if nan_counter >= 900 % 500
        
        if take_out_ind > 1
                    
                   add_row = cell(1,1);
                   take_out = vertcat(take_out, add_row);
                    
        end
        
        take_out{take_out_ind,1} = strcat('iEEG','_',num2str(elecs(e)));
        take_out_ind = take_out_ind + 1; 
             
    else
        
        if keep_in_ind > 1
                    
                   add_row = cell(1,1);
                   keep_in = vertcat(keep_in, add_row);
                    
        end
        
        keep_in{keep_in_ind,1} = strcat('iEEG','_',num2str(elecs(e)));
        keep_in_ind = keep_in_ind + 1; 

        
        ROL_Values = [ROL_Values, nanmedian(onsets(:,e))*1000];
        
    end
    
end

Post_Hoc_Elec_ROLs = cell(length(keep_in),2);

Post_Hoc_Elec_ROLs(:,1) = keep_in';
Post_Hoc_Elec_ROLs(:,2) = num2cell(round(ROL_Values)');

