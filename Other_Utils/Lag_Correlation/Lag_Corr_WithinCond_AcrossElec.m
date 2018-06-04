function [ATbT_Correlation_Values,lags] = Lag_Corr_WithinCond_AcrossElec(fname,indchan,conds,timewin, save,suffix,colors,labels)

% Function this function performes a lag correlation between two electrodes
% on the same condition
% Inputs:
% - fname     : name of file to display (SPM format)
% - indchan   : enter two channels
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

if nargin<1 || isempty(fname)
    fname = spm_select(1,'mat','Select file to display',{},pwd,'.mat');
end
D = spm_eeg_load(fname);

if nargin<2 || isempty(indchan)
    indchan = 1:D.nchannels;
elseif iscellstr(indchan)
    indchan = indchannel(D,indchan);
end

if nargin<3 || isempty(conds)
    labs = D.condlist;
elseif iscellstr(conds) || iscell(conds)
    labs = conds;
else
    listcond = condlist(D);
    labs  =listcond(conds);
end

if nargin<4 || isempty(timewin)
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

if nargin<5 || isempty(save)
    save = 1;
end

if nargin<6
    suffix = [];
end

if nargin<7  
   colors = colormap(lines(length(labs)-1));
   colors = [[200 10 150]/255; colors]; % pink
else
    if size(colors,1)< length(labs)
        disp('Number of colors provided smaller than number of conditions to plot')
        disp('Completing with Matlab built-in')
        addc = size(colors,1)-length(labs);
        colors = [colors; colormap(lines(addc))];
    end
end

if nargin<8
    labels = labs;
end

% check number of elecs
if length(labs) == 1 && length(indchan) == 2

    idt = tt(1):tt(2);

    tr_toplot = setdiff(indtrial(D,labs{1}),badtrials(D));

    chan1_mcond = squeeze(mean(D(indchan(1),idt,tr_toplot),3));
    chan2_mcond = squeeze(mean(D(indchan(2),idt,tr_toplot),3));

    % calculate across means
    [MXCF,lags,bounds] = crosscorr(chan1_mcond,chan2_mcond,500);
    bounds;

    ATbT_Correlation_Values = MXCF;


end











%%

% 
% chan1_concat = [];
% chan2_concat = [];

% 
% % calclulate across trials
% for i = 1:size(chan1_cond,3)
%     
%     [XCF,lags,bounds] = crosscorr(chan1_cond(1,:,i),chan2_cond(1,:,i),200);
%     bounds;
%     
%     % save concatinated info
%     chan1_concat = [chan1_concat, chan1_cond(1,:,i)];
%     chan2_concat = [chan2_concat, chan2_cond(1,:,i)];
%     
%     TbT_Correlation_Values(i,1:length(XCF)) = XCF;
%    
%     
% end



%  Plot the confidence bounds (horizontal lines) under the hypothesis that
%  the underlying series are uncorrelated.

% bin data 
% bins = length(MXCF)/101;
% 
% bin_MXCF = [];
% 
% 
% 
% for b = 1:101
%     
%     if b == 1
%         
%         eval(['bin_' num2str(b) '= MXCF(1:bins*b)']);
%         
%         
%     elseif b > 1 
%         
%         eval(['bin_' num2str(b) '= MXCF(bins*(b-1):bins*b)']);
%         
%     elseif b == bins
%         
%         eval(['bin_' num2str(b) '= MXCF(bins*(b-1):end)']);
%     
%     end
%     
%     eval(['bin_MXCF = [bin_MXCF, mean(bin_' num2str(b) ')];']);
%     
%     
% end
% 
% bin_lags = -1000:20:1000;
% 
% %% Plot
% 
% figure;
% 
% a = axis;
% 
% plot([1000 -1000],[0 0],'-k');
% 
% hold on
% 
%    lineHandles = bar(bin_lags,bin_MXCF,1,'FaceColor',[0,0.513725490196078,0.611764705882353],'EdgeColor',[0,0,0],'FaceAlpha',0.5);
%    %lineHandles(1).Color = ;
%    %set(lineHandles(1),'MarkerSize',3)
%    %set(lineHandles(1),'Color',[0.768627450980392,0.192156862745098,0.192156862745098])
%    %set(lineHandles(1),'LineWidth',1)
%    grid('off')
%    xlabel('Lag (ms)')
%    ylabel('Sample Cross Correlation')
%    %title('Sample Cross Correlation Function: Mean')
%    set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'LineWidth'   , 1         );
%    set(gca,'linewidth',1)
% 
%    
%    xlim([-200 1000])
%    ylim([-.6 0])
% 
%     a = get(gca,'XTickLabel');
%     set(gca,'XTickLabel',a,'fontsize',16,'fontweight','bold')
% 
%     y = get(gca,'XTickLabel');
%     set(gca,'XTickLabel',y,'fontsize',16,'fontweight','bold')
%     
%     [~,inx] = min(bin_MXCF);
%     m_lag =  bin_lags(inx);
%     hold on
%     plot([m_lag m_lag],[-6 0],'--k','LineWidth',2)
%     
%     
%     
% 
% ATbT_Correlation_Values = mean(TbT_Correlation_Values,1);
% 
%   
%    
% %  Plot the confidence bounds (horizontal lines) under the hypothesis that
% %  the underlying series are uncorrelated.
% 
% 
% % figure;
% % 
% % lineHandles = stem(lags,ATbT_Correlation_Values,'filled','r-o');
% %    set(lineHandles(1),'MarkerSize',4)
% %    grid('on')
% %    xlabel('Lag')
% %    ylabel('Sample Cross Correlation')
% %    title('Average across single trials')
% %    
% % figure;
% % 
% % 
% % [MXC,lags,bounds] = crosscorr(chan1_concat,chan2_concat,1000);
% % bounds;
% % 
% % lineHandles = stem(lags,MXC,'filled','r-o');
% %    set(lineHandles(1),'MarkerSize',4)
% %    grid('off')
% %    xlabel('Lag (ms)')
% %    ylabel('Sample Cross Correlation')
% %    title('over concat')
% 
% 
% 
% 
