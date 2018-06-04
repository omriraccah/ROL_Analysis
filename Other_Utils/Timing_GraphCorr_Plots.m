%% This script will be used to correlated trial-by-trial across regions and against RT

% Trial averages across electrodes in a specific region
PMC_Regional_Ave = nanmedian(tStruct.PMC.TbT_Timing, 1)*1000;
PFC_Regional_Ave = nanmedian(tStruct.PFC.TbT_Timing, 1)*1000;
AG_Regional_Ave = nanmedian(tStruct.AG.TbT_Timing, 1)*1000;

%% Plot Bar Graph

f = figure;

 
Y = [nanmean(PMC_Regional_Ave), ...
    nanmean(AG_Regional_Ave), ...
    nanmean(PFC_Regional_Ave)];
 
SEM1 = nanstd(PMC_Regional_Ave)/sqrt(length(PMC_Regional_Ave));
SEM2 = nanstd(AG_Regional_Ave)/sqrt(length(AG_Regional_Ave));
SEM3 = nanstd(PFC_Regional_Ave)/sqrt(length(PFC_Regional_Ave));

SEM = [SEM1 SEM2 SEM3];

C = [[0.317647058823529,0.607843137254902,0.662745098039216];
    [0,0.513725490196078,0.611764705882353];
    [0.105882352941176,0.376470588235294,0.592156862745098]];

P = nan(numel(Y), numel(Y));
% P(1,2) = 0.0165;
% P(5,4) = 1.99996000079998e-05;
% P(5,3) = 1.99996000079998e-05;
% P(5,2) = 1.99996000079998e-05;
% Make P symmetric, by copying the upper triangle onto the lower triangle
PT = P';
lidx = tril(true(size(P)), -1);
P(lidx) = PT(lidx);

hold on

b = superbar(Y,'BarWidth', .6   , 'E', SEM, 'P', P,'ErrorbarColor', [0 0 0],'BarFaceColor', C, ...
 'ErrorbarLineWidth', 2,'Orientation', 'v','PStarFontSize',25);

set(gcf,'Position',[500 500 350 450]);

ylabel('ROL (ms)','fontsize',16);
%ylim([0 .8])
title('S2','fontsize',16);

set(gca,'Fontsize',14,'FontWeight','bold','LineWidth',2,'TickDir','out');


 set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)
   
xticks([1 2 3])
set(gca,'XTickLabel',{'PMC', 'AG' ,'mPFC'},'FontSize',16);
%set(f, 'Position', [100, 100, 700, 300]);

%% Correlation Plots

%% PMC vs. PFC

r = figure 

plot(PMC_Regional_Ave,PFC_Regional_Ave,'LineStyle','none','Marker','o','MarkerFaceColor',[0.419607843137255,0.223529411764706,0.372549019607843] , 'MarkerEdgeColor', [0.419607843137255,0.223529411764706,0.372549019607843],'MarkerSize',10)

xlabel('PMC ROL (ms)','fontsize',16);
ylabel('mPFC ROL (ms)','fontsize',16);

set(gca,'Fontsize',14,'FontWeight','bold','LineWidth',2,'TickDir','out');

 set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)
   
hold on
l = lsline

[R,P] = corrcoef(PMC_Regional_Ave, PFC_Regional_Ave,'rows','complete');

title(strcat('R = ', num2str(R(2)),' P = ', num2str(P(2))))

set(l(1),'color',[0.419607843137255,0.223529411764706,0.372549019607843])

set(l(1),'LineWidth',2)

%% AG vs. PFC

r = figure 

plot(AG_Regional_Ave,PFC_Regional_Ave,'LineStyle','none','Marker','o','MarkerFaceColor',[0.419607843137255,0.223529411764706,0.372549019607843] , 'MarkerEdgeColor', [0.419607843137255,0.223529411764706,0.372549019607843],'MarkerSize',10)

xlabel('AG ROL (ms)','fontsize',16);
ylabel('mPFC ROL (ms)','fontsize',16);

set(gca,'Fontsize',14,'FontWeight','bold','LineWidth',2,'TickDir','out');

%xlim([100 700])
%ylim([100 1000])

 set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)
   
hold on
l = lsline

[R,P] = corrcoef(AG_Regional_Ave, PFC_Regional_Ave,'rows','complete');

title(strcat('R = ', num2str(R(2)),' P = ', num2str(P(2))))

set(l(1),'color',[0.419607843137255,0.223529411764706,0.372549019607843])

set(l(1),'LineWidth',2)

%% PMC vs. AG

r = figure 

plot(PMC_Regional_Ave,AG_Regional_Ave,'LineStyle','none','Marker','o','MarkerFaceColor',[0.419607843137255,0.223529411764706,0.372549019607843] , 'MarkerEdgeColor', [0.419607843137255,0.223529411764706,0.372549019607843],'MarkerSize',10)

xlabel('PMC ROL (ms)','fontsize',16);
ylabel('AG ROL (ms)','fontsize',16);

set(gca,'Fontsize',14,'FontWeight','bold','LineWidth',2,'TickDir','out');

%xlim([100 700])
%ylim([100 1000])

 set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)
   
hold on
l = lsline

[R,P] = corrcoef(PMC_Regional_Ave, AG_Regional_Ave,'rows','complete');

title(strcat('R = ', num2str(R(2)),' P = ', num2str(P(2))))

set(l(1),'color',[0.419607843137255,0.223529411764706,0.372549019607843])

set(l(1),'LineWidth',2)


