function Plot_two_freq_bands(fname1,fname2,indchan,conds,timewin, save,suffix,colors,~)

% Function to plot two different frequency bands (e.g. HFB and Theta) for a
% given condition on a single plot.
% Inputs:
% - fname     : name of file to display (SPM format)
% - indchan   : index or name of the channel to plot (default: 1st channel)
% - conds     : indexes of conditions to display or cell array of condition
%               names. Default: all conditions in file
% - timewin   : time window to display (in ms). Default: all time points
% - save      : save the plots in a 'figs' subdirectory of the file name 
%               path (Default: Yes, 1)
% suffix      : specific suffix for saving the plots
% colors      : color matrix for plot
% Outputs:
% Plots of the signal average and standard error for each condition. One
% plot per channel. Those plots can be saved automatically if wanted.
% -------------------------------------------------------------------------
% Written by J. Schrouff, LBCN, Stanford, 10/28/2015

% Get inputs
% -----------

if nargin<1 || isempty(fname1)
    fname = spm_select(1,'mat','Select file to display',{},pwd,'.mat');
end
D = spm_eeg_load(fname1);

if nargin<2 || isempty(fname2)
    fname2 = spm_select(1,'mat','Select file to display',{},pwd,'.mat');
end
D2 = spm_eeg_load(fname2);

if nargin<3 || isempty(indchan)
    indchan = 1:D.nchannels;
elseif iscellstr(indchan)
    indchan = indchannel(D,indchan);
end

if nargin<4 || isempty(conds)
    labs = D.condlist;
elseif iscellstr(conds) || iscell(conds)
    labs = conds;
else
    listcond = condlist(D);
    labs  =listcond(conds);
end

if nargin<5 || isempty(timewin)
    time_start = min(time(D));
    tt(1) = indsample(D,time_start);
    time_end = max(time(D));
    tt(2) = indsample(D,time_end);
else
    tt(1) = indsample(D,timewin(1)/1000);
    tt(2) = indsample(D,timewin(2)/1000);
    time_start = time(D,tt(1));
    time_end = time(D,tt(2));
end

if nargin<6 || isempty(save)
    save = 1;
end

if nargin<7
    suffix = [];
end

if nargin<8  
   colors = colormap(lines(length(labs)));
   colors = [[200 10 150]/255; colors]; % pink
else
    if size(colors,1)< length(labs)
        disp('Number of colors provided smaller than number of conditions to plot')
        disp('Completing with Matlab built-in')
        addc = size(colors,1)-length(labs);
        colors = [colors; colormap(lines(addc))];
    end
end

if nargin<9
    labels = labs;
end


% Plot traces
% -----------

tim =time(D);
idt = tt(1):tt(2);

if numel(size(D)) == 4
    ifreq = 1;
    disp('Plotting for first frequency')
    itf = 1;
else
    itf = 0;
end

plot_max = zeros(length(labs),1);
plot_min = zeros(length(labs),1);

fprintf(['Plotting channel (out of %d):',repmat(' ',1,ceil(log10(length(indchan)))),'%d'],length(indchan), 1);
for i=1:length(indchan)
    
    % Counter of channels to be updated
    if i>1
        for idisp = 1:ceil(log10(i)) % delete previous counter display
            fprintf('\b');
        end
        fprintf('%d',i);
    end
            
    hfig =figure;
    set(gcf,'Position',[400 500 800 500]);
    hold on
    
    % First plot the signal to create the legend
    for sp = 1:2
        tr_toplot = setdiff(indtrial(D,labs{1}),badtrials(D)); %take bad trials out
        if itf
            mcond = squeeze(mean(D(indchan(i),itf,idt,tr_toplot),4));
        else
            mcond = squeeze(mean(D(indchan(i),idt,tr_toplot),3));
        end
        plot(tim(idt),mcond, 'Linewidth',2,'Color',colors(sp,:));
    end    
    legend({'Theta','HFB'},'Location','NorthEastOutside'); 
    
    % Add standard error for each channel
    for sp = 1:length(labs)   
        tr_toplot = setdiff(indtrial(D,labs{sp}),badtrials(D));
         if itf
            mcond = squeeze(mean(D(indchan(i),itf,idt,tr_toplot),4));
            stdcond = squeeze(std(D(indchan(i),itf,idt,tr_toplot),1,4))/sqrt(length(tr_toplot));
        else
            mcond = squeeze(mean(D(indchan(i),idt,tr_toplot),3));
            stdcond = squeeze(std(D(indchan(i),idt,tr_toplot),1,3))/sqrt(length(tr_toplot));
        end

        if ~isempty(~isnan(mcond))
            shadedErrorBar(tim(idt),mcond,stdcond,{'linewidth',2,'Color',colors(1,:)},0.8);
            %for axis below
            plot_max(sp,:) = max(mcond+stdcond);
            plot_min(sp,:) = min(mcond-stdcond);
        end
        hold on
        
         tr_toplot2 = setdiff(indtrial(D2,labs{sp}),badtrials(D2));
         if itf
            mcond2 = squeeze(mean(D2(indchan(i),itf,idt,tr_toplot2),4));
            stdcond2 = squeeze(std(D2(indchan(i),itf,idt,tr_toplot2),1,4))/sqrt(length(tr_toplot2));
        else
            mcond2 = squeeze(mean(D2(indchan(i),idt,tr_toplot2),3));
            stdcond2 = squeeze(std(D2(indchan(i),idt,tr_toplot2),1,3))/sqrt(length(tr_toplot2));
        end

        if ~isempty(~isnan(mcond))
            shadedErrorBar(tim(idt),mcond2,stdcond2,{'linewidth',2,'Color',colors(2,:)},0.8);
            %for axis below
            plot_max(sp,:) = max(mcond2+stdcond2);
            plot_min(sp,:) = min(mcond2-stdcond2);
        end
        
    end
    
    
    %add lines
    %zero level
    line([time_start time_end],[0 0],'LineWidth',1,'Color','k');
    % %stim onset
    %line([0 0],[min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10],'LineWidth',1,'Color','k'); %stim onset;
    
    
    %figure properties
    hold off
    box off
    xlim([time_start time_end]);
    %ylim([min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10]);
    xlabel('Time (s)','fontsize',14);
    ylabel('Signal','fontsize',14);
    if ismember(indchan(i),D.badchannels)
        isgood = 'bad';
    else
        isgood = 'good';
    end
    title(['Channel ',char(chanlabels(D,indchan(i))),'  (',isgood,')'],'fontsize',14, 'fontweight','bold')
    
    % Save figures in a 'figs' subdirectory
    path = pwd;
    if save
        if ~exist([path,filesep,'figs'],'dir')
            mkdir(path,'figs');
        end
        print('-opengl','-r300','-dpng',strcat([path,filesep,'figs',filesep,char(chanlabels(D,indchan(i))),'_',suffix]));
        delete(hfig)
    end
    
end
fprintf('\n')
disp('Done: Plot average and standard error for each channel')