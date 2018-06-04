% Prepare files and plot spectrogram for gradCPT RT variability

% Written by Aaron Kucyi and Omri Raccah, LBCN, Stanford University
%==========================================================================

% %% Epoch mountains
% display(['Choose preprocessed/re-referenced file & Mountain events file']);
% eMount=LBCN_epoch_bc([],[],[],'start',[-2500 2500],1,'start',[-800 -50],'eMount');

function [] = Spectograms(spm_path)



%% Convert from SPM to field trip
eMount_ft=spm2fieldtrip(spm_file);

%% Time-frequency decomposition in fieldtrip
cfg=[];
cfg.output='pow';
cfg.method='wavelet';
freq = [1 200];
     nf=length(freq(1):3:freq(2));
     cfg.foi=logspace(log10(freq(1)),log10(freq(2)),nf);
cfg.toi=-0.9:0.05:0.9;
cfg.keeptrials='yes';
tf_eMount_ft=ft_freqanalysis(cfg,eMount_ft);
save('tf_eMount_ft','tf_eMount_ft');

%% Baseline normalize
cfg.baseline=[-0.8 -0.05];
cfg.baselinetype='relative';
btf_eMount_ft=ft_freqbaseline(cfg,tf_eMount_ft);

%% Average across trials
cfg=[];
cfg.parameter='powspctrm';
avg_btf_eMount_ft=ft_freqgrandaverage(cfg,tf_eMount_ft);

%% Display channel numbering/naming in matlab window
chan_nums=num2cell((1:length(btf_eMount_ft.label))');
chan_labels=btf_eMount_ft.label;
[chan_nums chan_labels]

%% Choose channel to plot
Channels=input('Channels to plot, e.g. [51 52]: ');

%% Plot spectrograms and save to file
mkdir(['figs/Spectrograms']);
cd(['figs/Spectrograms']);

for i=1:length(Channels)
    chanlabel=char(avg_btf_eMount_ft.label(Channels(i)));
    
figure(1);
contourf(avg_btf_eMount_ft.time, avg_btf_eMount_ft.freq,squeeze(avg_btf_eMount_ft.powspctrm(Channels(i),:,:)),40,'linestyle','none');
axis([-0.9 0.9 2 200]);
colormap('hot'); colorbar;
set(gca,'yscale','log','ytick',[2 4 8 16 32 64 128]);
%caxis([0.5 3]);
xlabel('Time (sec)'); ylabel('Frequency (Hz)');
set(gcf,'PaperPositionMode','auto');
print([chanlabel '_Mountains'],'-r300','-dpng');

end
close('all');