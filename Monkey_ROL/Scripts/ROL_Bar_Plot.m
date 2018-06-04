function ROL_Bar_Plot(fname,indchan,cond,ROL,bin_data,ROL_Bin_Index,sig_thr,timewin,save,prefix,colors)

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

set(0,'DefaultFigureVisible','off')

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

if nargin<6 || isempty(timewin)
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
if nargin<9 || isempty(save)
    save = 1;
end

if nargin<10
    prefix = [];
end

if nargin<11
   colors = colormap(lines(5));
else
    if size(colors,1)< length(indchan)
        disp('Number of colors provided smaller than number of channels to plot')
        disp('Completing with Matlab built-in')
        addc = size(colors,1)-length(indchan);
        colors = [colors; colormap(lines(addc))];
    end
end

%% Plot bar
% -----------
% re-create figure lables    
c_names = {};
for cn = 1:length(indchan)
    
    figure
    
    % x-values (onset_indeces)
    time_values =  bin_data{1,cn}(:,5);
    time_values = cell2mat(time_values);
    time_values = D.time(time_values)*1000;
    ROL_time = time_values(ROL_Bin_Index{1, 1});
    
    onset_str =  [' ,onset: ', num2str(ROL_time)]; 
    chan_str =  ['channel: ', num2str(indchan(cn))];   
    
    name =  strcat(chan_str, onset_str);

    c_names = [c_names, name];

end

for i= 1:length(indchan)
    
    f = figure;
    
    % x-values (onset_indeces)
    time_values =  bin_data{1,i}(:,5);
    time_values = cell2mat(time_values);
    time_values = D.time(time_values)*1000;
    ROL_time = time_values(ROL_Bin_Index{1, i});
    
    % y-values (means)
    mean_bin_values = bin_data{1,i}(:,2);
    mean_bin_values = cell2mat(mean_bin_values);
    ROL_mean = mean_bin_values(ROL_Bin_Index{1, i});
    
    slopes = cell2mat(bin_data{1,i}(:,3));
    ROL_slope = slopes(ROL_Bin_Index{1, i})
    
    pos_slope_indx = find(slopes >= 0);
    neg_slope_indx = find(slopes < 0);
    
    p1 = bar(time_values(pos_slope_indx),mean_bin_values(pos_slope_indx));
    hold on;
    p2 =  bar(time_values(neg_slope_indx),mean_bin_values(neg_slope_indx));
    set(p1,'FaceColor', colors(2,:));
    set(p2,'FaceColor',colors(1,:));
    plot(xlim,[sig_thr{1,i} sig_thr{1,i}],'--','LineWidth',2,'Color','k') 
   
    hold on
    % mark ROL
    if ROL_mean > 0
        plot(ROL_time,ROL_mean+.05,'*','MarkerSize', 12,'Color', [0 0.7  0],'linewidth',2)
        
    else 
        plot(ROL_time,ROL_mean-.05,'*','MarkerSize', 12,'Color', [0 0.7  0],'linewidth',2)
    
    end
    
    xbound = xlim;
    
    xlim([0,800])
    
    xlabel('Time (ms)','fontsize',20);
    ylabel('Signal','fontsize',20);
    set(f, 'Position', [1000, 1000, 4049, 1395]);

    % title(c_names{1, i},'fontsize',35, 'fontweight','bold')

    % Save figures in a 'figs' subdirectory in current working path
    if save
       if ~exist([pwd,filesep,'figs'],'dir')
                mkdir(pwd,'figs');
       end
            print('-dpng',strcat([pwd,filesep,'figs',filesep,prefix,'Monkey_ROL_MeanSignal_', fname(end-11: end-4),'_channel_',num2str(indchan(i))]));
    end

fprintf('\n')
disp('Done: Plot average and standard error for each channel')

    
end




