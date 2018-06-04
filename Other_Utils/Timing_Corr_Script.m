%% This script will be used to correlated trial-by-trial across regions and against RT

% Trial averages across electrodes in a specific region
PMC_Regional_Ave = nanmedian(tStruct.PMC.TbT_Timing, 1);
PFC_Regional_Ave = nanmedian(tStruct.PMC.TbT_Timing, 1);
AG_Regional_Ave = nanmedian(tStruct.AG.TbT_Timing, 1);




