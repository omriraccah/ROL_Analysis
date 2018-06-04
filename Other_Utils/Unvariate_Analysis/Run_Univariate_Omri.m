addpath /Users/parvizilab/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_38_LK/MMR_SPM_Files/LK_08

fname1 = '/Users/parvizilab/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_38_LK/MMR_SPM_Files/LK_08/SHFBrtf_aeMfffspm8_iEEGLK_08.mat';

LK_baseline = LBCN_Univariate_Condition_Baseline_OR(fname1,{'math'},[100 1000],[-100 0],'good',50000,[]);

LK_condition = LBCN_Univariate_Condition_Compare_OR(fname1,{'math'},[400 800],{'autobio'},[400 800],'good',10000,[],[]);


fname2 = '/Users/labuser/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_42_NC/MMR_SPM_Files/NC_06/SHFBrtf_aeMfffspm8_iEEGNC_06.mat';

NC_baseline = LBCN_Univariate_Condition_Baseline_OR(fname2,{'math'},[100 1000],[-100 0],'good',10000,[]);

NC_condition = LBCN_Univariate_Condition_Compare_OR(fname2,{'math'},[400 800],{'autobio'},[400 800],'good',10000,[],[]);