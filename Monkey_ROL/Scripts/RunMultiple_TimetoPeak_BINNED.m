%% Data paths

paths = {'/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Preprocessed_Data/S12_38_LK/StimulusBased/cHFBrtf_aeMfffspm8_iEEGLK_08.mat'...
    ,'/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Preprocessed_Data/S12_42_NC/Stimulus_Based/cHFBrtf_aeMfffspm8_iEEGNC_05.mat'...
    ,'/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Preprocessed_Data/S13_47_JT2/Stimulus_Based/cHFBrtf_aeMfffspm8_iEEGJT2_01.mat'...
    ,'/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Preprocessed_Data/S17_110_SC/Stimulus_Based/cHFBrtf_aeMfffECoG_E17-394_0006.mat'...
    ,'/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Preprocessed_Data/S17_111_RT/Stimulus_Based/cHFBrtf_aeMfffECoG_E17-438_0008.mat'...
    ,'/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Preprocessed_Data/S11_28_SRb/cHFBrtf_aeMfffspm8_iEEGSRb-06.mat'};

%% Math

S1_Math = {'iEEG_66','iEEG_67','iEEG_82','iEEG_83','iEEG_84','iEEG_89','iEEG_43','iEEG_52','iEEG_59','iEEG_44','iEEG_60','iEEG_34','iEEG_42','iEEG_49'}; S1_deactive = [1,1,1,1,1,1,0,0,0,0,0,0,0,0];
S2_Math = {'iEEG_73','iEEG_87','iEEG_39','iEEG_38','iEEG_40','iEEG_55','iEEG_33','iEEG_34','iEEG_35'}; S2_deactive = [1,1,0,0,0,0,0,0,0];
S3_Math = {'iEEG_55','iEEG_57','iEEG_107','iEEG_44'}; S3_deactive = [1,1,0,0];
S4_Math = {'LPC1','LPC7','LPC8','RO2','RO3','RO4','RO5','RO6','RO8'}; S4_deactive = [1,0,0,0,0,0,0,0,0];
S5_Math = {'LP2','LP3','LPI12','LPS1','LPI17','LPI19','LPI48','LP7','LP8','LP10','LP6'}; S5_deactive = [1,1,1,1,0,0,0,0,0,0,0];
S6_Math = {'iEEG_74','iEEG_12','iEEG_20'}; S6_deactive = [1,0,0];

allElecs = {S1_Math;S2_Math;S3_Math;S4_Math;S5_Math;S6_Math};
allDeactive = {S1_deactive;S2_deactive;S3_deactive;S4_deactive;S5_deactive;S6_deactive};

Peak_Results = struct;

for s = 1:length(paths)

         deactive_arr = allDeactive{s};        
        
         peak = TbT_Main(paths{s},allElecs{s, 1},'math',1, allDeactive{s},0);

            
         eval(['Peak_Results.S' num2str(s) '.math{c,1} = elec;'])
         eval(['Peak_Results.S' num2str(s) '.math{c,2} = round(peak*1000);'])
        
        
        disp(strcat('done with subject:', num2str(s)))

end
    

%% Autobio

% 
% S1_Auto = {'iEEG_66','iEEG_67','iEEG_68','iEEG_71','iEEG_72','iEEG_82','iEEG_83','iEEG_84' ,'iEEG_86','iEEG_87','iEEG_88','iEEG_89','iEEG_34','iEEG_42','iEEG_49','iEEG_59','iEEG_44','iEEG_60'}; S1_deactive = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% S2_Auto = {'iEEG_73','iEEG_74','iEEG_83','iEEG_84','iEEG_38','iEEG_40','iEEG_55','iEEG_33','iEEG_34','iEEG_35'}; S2_deactive = [0,0,0,0,0,0,0,0,0,0];
% S3_Auto = {'iEEG_54','iEEG_57','iEEG_58','iEEG_44'}; S3_deactive = [0,0,0,0];
% S4_Auto = {'RO2','RO3','RO4','RO5','RO6','RO8'}; S4_deactive = [0,0,0,0,0,0];
% S5_Auto = {'LP1','LP2','LP4','LPI11','LPI12','LPS1','LPI48','LP7','LP8','LP10','LP6'}; S5_deactive = [0,0,0,0,0,0,0,0,0,0,0];
% S6_Auto = {'iEEG_74','iEEG_20'}; S6_deactive = [0,0];
% 
% allElecs = {S1_Auto;S2_Auto;S3_Auto;S4_Auto;S5_Auto;S6_Auto};
% allDeactive = {S1_deactive;S2_deactive;S3_deactive;S4_deactive;S5_deactive;S6_deactive};
% 
% for s = 1:length(paths)
% 
%      deactive_arr = allDeactive{s};
%     
%     % loop through electrodes
%     for c = 1:length(allElecs{s,1})
%         
%         elec = allElecs{s, 1}{1, c};
%         D = spm_eeg_load(paths{s});
%         
%         if s == 4 || s == 5
%            
%          elec = indchannel(D,elec);
%           
%         else 
%            
%          elec = str2num(elec(6:end));
%         
%         end
%         
%         indtr =  setdiff(indtrial(D,'autobio'),badtrials(D));
%         
%         
%          peak = Fit_T2P_WideBins(D,elec,indtr,[], deactive_arr(c),0);
% 
%             
%          eval(['Peak_Results.S' num2str(s) '.auto{c,1} = elec;'])
%          eval(['Peak_Results.S' num2str(s) '.auto{c,2} = round(peak*1000);'])
%             
%             
%     end
%         
%         
%         disp(strcat('done with subject:', num2str(s)))
% 
% end