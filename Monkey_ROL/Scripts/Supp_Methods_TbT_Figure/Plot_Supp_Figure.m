function [] = Plot_Supp_Figure(path,elec,trial,condition,range)

%% Load trial data

% load subject data
D = spm_eeg_load(path);

tt(1) = indsample(D,range(1)/1000);
tt(2) = indsample(D,range(2)/1000);
idt = tt(1):tt(2);

% Get all good trials
% tr_toplot = setdiff(indtrial(D,condition),badtrials(D));

% Create matrix (trials vs. time points)
data = squeeze(D(elec,idt,trial));

%% Plot trial data

plot(D.time(idt), data,'LineWidth',2)


   
    %add lines
    %zero level
    line([time_start time_end],[0 0],'LineWidth',2,'Color','k');
    % %stim onset
    ylim([-.3 1.6])
    line([0 0],ylim,'LineWidth',2,'Color','k'); %stim onset;
    
    
    %figure properties
    hold off
    box off
    xlim([-200 1000]);
    line(xlim,[0.0 0],'LineWidth',2,'Color','k');
    ylabel('HFB (z-score)','fontsize',16);

    
    set(gca,'Fontsize',16,'FontWeight','bold','LineWidth',2,'TickDir','out');


    set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)
 set(gca, 'FontName', 'Arial')

    
    %ylim([min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10]);
    xlabel('Time (ms)','fontsize',16);
    ylabel('HFB (z-score)','fontsize',16);