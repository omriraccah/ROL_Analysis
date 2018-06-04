function Plot_Common_Elec_Normalized_Signals(fname,indchan,cond,ROL,timewin,save,prefix,colors,deactive)

% Function to plot the signals from multiple electrodes (from a single
% subject) on a common plot for a specific condition (e.g. faces)
% Inputs:
% - fname     : name of file to display (SPM format)
% - indchan   : index or name of the channel to plot (default: 1st channel)
% - cond     :  condition of interest (default: 1st condition)
% - ROL     :  this is the ROL matrix gathered in the analysis scripts (no default)
%    
% - timewin   : time window to display (in ms). Default: all time points
% - save      : save the plots in a 'figs' subdirectory of the file name 
%               path (Default: Yes, 1)
% prefix      : specific prefix for saving the plots
% colors      : color matrix for plot
% Outputs:
% Plots of the signal average for all electrodes on a signal plot.Those plots can be 
% saved automatically unless specified otherwise.
% -------------------------------------------------------------------------
% Written by O. Raccah (adapted from J. Schrouff's 'LBCN_plot_average_signal_epochs' preprocessing script), LBCN, Stanford, 09/08/2016

% Get inputs
% -----------


if nargin<1 || isempty(fname)
    fname = spm_select(1,'mat','Select file to display',{},pwd,'.mat');
end
D = spm_eeg_load(fname);

if nargin<2 || isempty(indchan)
    indchan = 1:D.nchannels;
elseif iscellstr(indchan)
    indchan = indchannel(D,indchan);
end

if nargin<3 || isempty(cond)
    labs = D.condlist(1);
elseif iscellstr(cond)
    labs = cond;
else
    listcond = condlist(D);
    labs  =listcond(cond);
end

if nargin<5 || isempty(timewin)
    tt(1) = indsample(D,-.1);
    tt(2) = indsample(D,1.8);
    time_start = time(D,tt(1));
    time_end = time(D,tt(2));
else
    tt(1) = indsample(D,timewin(1)/1000);
    tt(2) = indsample(D,timewin(2)/1000);
    time_start = time(D,tt(1));
    time_end = time(D,tt(2));
end

if nargin<6 || isempty(save)
    save = 1;
end

if nargin<7|| isempty(save)
    prefix = [];
end

if nargin<8|| isempty(colors)
   colors = colormap(lines(length(indchan)));
else
    if size(colors,1)< length(indchan)
        disp('Number of colors provided smaller than number of channels to plot')
        disp('Completing with Matlab built-in')
        addc = size(colors,1)-length(indchan);
        colors = [colors; colormap(lines(addc))];
    end
end

%% Plot traces
% -----------

set(0,'DefaultFigureVisible','off')

tim =time(D)*1000;
idt = tt(1):tt(2);

plot_max = zeros(length(labs),1);
plot_min = zeros(length(labs),1);

p_all = [];

% save channel names

c_names = {};
for cn = 1:length(indchan)
    
    c_names = [c_names, num2str(indchan(cn))];

end
 

hfig =figure;
set(gcf,'Position',[400 500 800 500]);
    

for i= 1:length(indchan)
    
    % Counter of channels to be updated
    if i>1
        for idisp = 1:ceil(log10(i)) % delete previous counter display
            fprintf('\b');
        end
        fprintf('%d',i);
    end
    
    hold on
           
    % First plot the signal to create the legend
        tr_toplot = setdiff(indtrial(D,labs),badtrials(D)); %take bad trials out
        
        % ************ Normalize by Trial
        
        norm_data = D(indchan(i),idt,tr_toplot);
    
        if deactive(i) == 0
            for n = 1:length(norm_data(1,1,:))

                MT = max(norm_data(:,:,n));

                norm_data(:,:,n) = norm_data(:,:,n) / MT;

            end

        elseif deactive(i) == 1
            
            for n = 1:length(norm_data(1,1,:))

                MT = min(norm_data(:,:,n));

                norm_data(:,:,n) = norm_data(:,:,n) / abs(MT);

            end
            
            
        end
        
        
        % ************
        
        
        
        mcond = mean(norm_data,3);
        
        eval(['p_' num2str(i) '=' 'plot(tim(idt),mcond, ''Linewidth'',2,''Color'',colors(i,:));'])
        eval(['p_all = [p_all, p_' num2str(i) '];'])

%     % Add standard error for each channel
%         tr_toplot = setdiff(indtrial(D,labs),badtrials(D));
%         mcond = mean(D(indchan(i),idt,tr_toplot),3);
        
       % stdcond = std(mcond,1,3)/sqrt(length(tr_toplot));
        
        stdcond = std(norm_data,1,3)/sqrt(length(tr_toplot)); %
        %commented for normalization ^^
        if ~isempty(~isnan(mcond))
            shadedErrorBar(tim(idt),mcond,stdcond,{'linewidth',2,'Color',colors(i,:)},0.8);
            %for axis below
            plot_max(i,:) = max(mcond+stdcond);
            plot_min(i,:) = min(mcond-stdcond);
        end
        hold on
    
    %add lines
    %zero level
    line([time_start time_end],[0 0],'LineWidth',1,'Color','k');
    % %stim onset
    line([0 0],[min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10],'LineWidth',1,'Color','k'); %stim onset;
    
    
    %figure properties
    hold off
    box off
    xlim([-100 1800]);
    ylim([min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10]);
    xlabel('Time (ms)','fontsize',14);
    ylabel('Signal','fontsize',14);
    if ismember(indchan(i),D.badchannels)
        isgood = 'bad';
    else
        isgood = 'good';
    end
    title('All Face Selective Electrodes','fontsize',14, 'fontweight','bold')

end

legend(p_all,c_names,'Location','NorthEastOutside');


% Save figures in a 'figs' subdirectory in current working path
if save
   if ~exist([pwd,filesep,'figs'],'dir')
            mkdir(pwd,'figs');
   end
        print('-opengl','-r300','-dpng',strcat([pwd,filesep,'figs',filesep,prefix,'SignalNormalized_ALL_conditions_', fname(end-11: end-4)]));
        delete(hfig)
end

fprintf('\n')
disp('Done: Plot average and standard error for each channel')


%% Plot with onsets without line error bars 
% -----------

tim =time(D)*1000;
idt = tt(1):tt(2);

nfig =figure;
set(gcf,'Position',[400 500 800 500]);
p_all = [];

% c_names = {};
% for cn = 1:length(indchan)
%     
%     c_names = [c_names, num2str(indchan(cn))];
% 
% end
 

% re-create figure lables    
c_names = {};
for cn = 1:length(indchan)
    
    onset_str =  [', onset: ', num2str(nanmedian(ROL.onsets(:,cn)*1000))]; 
    chan_str =  ['channel: ', num2str(indchan(cn))];   
    
    name =  strcat(chan_str, onset_str);

    c_names = [c_names, name];

end

%fprintf(['Plotting channel (out of %d):',repmat(' ',1,ceil(log10(length(indchan)))),'%d'],length(indchan), 1);
for ii= 1:length(indchan)
    
    % Counter of channels to be updated
    if ii>1
        for idisp = 1:ceil(log10(ii)) % delete previous counter display
            fprintf('\b');
        end
        fprintf('%d',ii);
    end
    
    hold on
           
    % First plot the signal to create the legend
        tr_toplot = setdiff(indtrial(D,labs),badtrials(D)); %take bad trials out
        
        % ************ Normalize by Trial
        
        norm_data = D(indchan(ii),idt,tr_toplot);
    
        if deactive(i) == 0
            for n = 1:length(norm_data(1,1,:))

                MT = max(norm_data(:,:,n));

                norm_data(:,:,n) = norm_data(:,:,n) / MT;

            end

        elseif deactive(i) == 1
            
            for n = 1:length(norm_data(1,1,:))

                MT = min(norm_data(:,:,n));

                norm_data(:,:,n) = norm_data(:,:,n) / abs(MT);

            end
            
            
        end
        
        mcond = mean(norm_data,3);
        
        
        eval(['p_' num2str(i) '=' 'plot(tim(idt),mcond, ''Linewidth'',2,''Color'',colors(ii,:));'])
        eval(['p_all = [p_all, p_' num2str(i) '];'])

        
    %add lines
    %zero level
    line([time_start time_end],[0 0],'LineWidth',1,'Color','k');
    % %stim onset
    line([0 0],[min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10],'LineWidth',1,'Color','k'); %stim onset;
    
    % ***** Plot ROL onsets for each electrode ******** 
    
    % calulate 95% confidence interval
    SE= nanmedian(ROL.onsets(:,ii)*1000)/sqrt(1000); %% Standard Error
    ME= (SE) *2 ;       %% margine of Error
    CI= ([ME+nanmedian(ROL.onsets(:,ii)*1000) , nanmedian(ROL.onsets(:,ii)*1000)-ME]);

    % % signal plot

    hold on
    plot([nanmedian(ROL.onsets(:,ii)*1000) nanmedian(ROL.onsets(:,ii)*1000)],ylim, 'Color', colors(ii,:),'LineWidth',1.2) 

    
    y = ylim;
    
    % Add a patch
    patch([CI(2), CI(1), CI(1), CI(2)],[y(1), y(1), y(2), y(2)],colors(ii,:),'EdgeColor',colors(ii,:),'EdgeAlpha','0.1')
    alpha(0.1)

    % ********* end **********
    
    
    %figure properties
    hold off
    box off
    xlim([-100 1800]);
    ylim([min(plot_min)-min(plot_min)*-.10 max(plot_max)+max(plot_max)*.10]);
    xlabel('Time (ms)','fontsize',14);
    ylabel('Signal','fontsize',14);
    if ismember(indchan(ii),D.badchannels)
        isgood = 'bad';
    else
        isgood = 'good';
    end
    title('Onsets for Selective Electrodes','fontsize',14, 'fontweight','bold')
    
end

% Plot legend with time values 
legend(p_all,c_names,'Location','NorthEastOutside');


% Save figures in a 'figs' subdirectory
%path = spm_fileparts(fname);
if save
   if ~exist([pwd,filesep,'figs'],'dir')
            mkdir(pwd,'figs');
   end
        print('-opengl','-r300','-dpng',strcat([pwd,filesep,'figs',filesep,prefix,'Signal_Normalized_ALLWithOnsets_conditions_', fname(end-11: end-4)]));
        delete(nfig)
end

fprintf('\n')
disp('Done: Plot average and standard error for each channel')



