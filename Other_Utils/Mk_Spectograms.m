% Prepare files and plot spectrogram for gradCPT RT variability

% Written by Aaron Kucyi and Omri Raccah, LBCN, Stanford University
%==========================================================================

% %% Epoch mountains
% display(['Choose preprocessed/re-referenced file & Mountain events file']);
% eMount=LBCN_epoch_bc([],[],[],'start',[-2500 2500],1,'start',[-800 -50],'eMount');

% function [] = Mk_Spectograms(spm_path)

%S1

%eMount=LBCN_epoch_bc('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_38_LK/MMR_TDT_Files/THIS_IS_LK_FInal/Mfffspm8_iEEGLK_10.mat','/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_38_LK/MMR_TDT_Files/THIS_IS_LK_FInal/eventsSODATA_spm8_iEEGLK_10.mat',5,'start',[-20000 4000],0,'start',[-200 0],'math');
% Merge_File = spm_eeg_load('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_38_LK/MMR_TDT_Files/THIS_IS_LK_FInal/cmathMfffspm8_iEEGLK_08.mat');


% S2 
%eMount=LBCN_epoch_bc('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_42_NC/MMR_SPM_Files/NC_Combined/Mfffspm8_iEEGNC_06.mat','/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_42_NC/MMR_SPM_Files/NC_Combined/eventsSODATA_spm8_iEEGNC_06.mat',4,'start',[-20000 4000],0,'start',[-200 0],'autobio');
%Merge_File = spm_eeg_load('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_42_NC/MMR_SPM_Files/NC_Combined/cautobioMfffspm8_iEEGNC_05.mat');

% % S3
%eMount=LBCN_epoch_bc('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S17_111_RT/MMR/MfffECoG_E17-438_0010.mat','/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S17_111_RT/MMR/eventsSODATA_DCchans_E17-438_0010.mat',4,'start',[-20000 4000],0,'start',[-200 0],'autobio');
% Merge_File = spm_eeg_load('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S12_38_LK/MMR_TDT_Files/THIS_IS_LK_FInal/cautobioMfffECoG_E17-438_0008.mat');


% % S4
% eMount=LBCN_epoch_bc('/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S13_47_JT2/MMR_SPM_RTBased/Mfffspm8_iEEGJT2_01.mat','/Users/oraccah/Documents/Omri/Network_Dynamics_Project/Unprocessed_Data/S13_47_JT2/MMR_SPM_RTBased/eventsSODATA_spm8_iEEGJT2_01.mat',4,'start',[-4000 4000],0,'start',[-200 0],'autobio');
% 

%% Convert from SPM to field trip
eMount_ft=spm2fieldtrip(Merge_File);

%% Time-frequency decomposition in fieldtrip
cfg=[];
cfg.output='pow';
cfg.method='wavelet';
freq = [1 200];
     nf=length(freq(1):3:freq(2));
     cfg.foi=logspace(log10(freq(1)),log10(freq(2)),nf);
cfg.toi= -.2:.05:2;
%cfg.keeptrials='yes';
tf_eMount_ft=ft_freqanalysis(cfg,eMount_ft);
save('tf_eMount_ft','tf_eMount_ft');

%% Baseline normalize
cfg.baseline=[-20 0];
cfg.baselinetype='db';
btf_eMount_ft=ft_freqbaseline(cfg,tf_eMount_ft);

%% Average across trials
cfg=[];
cfg.parameter='powspctrm';
avg_btf_eMount_ft=ft_freqgrandaverage(cfg,btf_eMount_ft);

%% Display channel numbering/naming in matlab window
chan_nums=num2cell((1:length(btf_eMount_ft.label))');
chan_labels=btf_eMount_ft.label;
[chan_nums chan_labels]

%% Choose channel to plot
Channels=input('Channels to plot, e.g. [51 52]: ');

%% Plot spectrograms and save to file
% mkdir(['figs/Spectrograms']);
% cd(['figs/Spectrograms']);

figure
%subplot(2,5,5)

for i=1:length(Channels)
    chanlabel=char(avg_btf_eMount_ft.label(Channels(i)));
    
    
%figure    
subplot(1,1,i)
contourf(avg_btf_eMount_ft.time, avg_btf_eMount_ft.freq,squeeze(avg_btf_eMount_ft.powspctrm(Channels(i),:,:)),40,'linestyle','none');
axis([-.2 1.5 2 200]);
colormap(cm); 
% 
% if i == 3
% colorbar;
% 
% end

set(gca,'yscale','log','ytick',[2 4 8 16 32 64 128]);
caxis([-3 3]);
shading 'interp'
% xlabel('Time (sec)','fontsize',16); %ylabel('Frequency (Hz)','fontsize',10);
%set(gcf,'PaperPositionMode','auto');

 set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',1)
   
hold on
plot([0 0],ylim,'-k','linewidth',1)


%print([chanlabel '_S1_autobio_chan_nocolorbar_sub' num2str(Channels(i))],'-r300','-dpng');

end
%close('all');