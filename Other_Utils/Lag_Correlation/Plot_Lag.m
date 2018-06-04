within = vertcat(within_66,within_67,within_82,within_83,within_84,within_89);

withinm = mean(within,1);

figure;

a = axis;

plot([1500 -1500],[0 0],'-k');


figure;

a = axis;

plot([1500 -1500],[0 0],'-k');

hold on

   lineHandles = stem(lags,within_82);
   %lineHandles(1).Color = ;
   set(lineHandles(1),'MarkerSize',5)
   set(lineHandles(1),'Color',[0.419607843137255,0.384313725490196,0.603921568627451])
   set(lineHandles(1),'LineWidth',2)
   grid('off')
   xlabel('Lag')
   ylabel('Sample Cross Correlation')
   %title('PMC: Autobio vs. Math (Across Electrodes)')
   set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)

   
   

    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',18,'fontweight','bold')

    y = get(gca,'XTickLabel');
    set(gca,'XTickLabel',y,'fontsize',18,'fontweight','bold')

