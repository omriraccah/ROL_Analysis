function [] = Plot_RT_HeatMap(path,elec,condition,range,deactive,interp,cm_new)

%% Load subject data and create trial vs. data matrix

% load subject data
D = spm_eeg_load(path);

behav = load(path);

tt(1) = indsample(D,range(1)/1000);
tt(2) = indsample(D,range(2)/1000);
idt = tt(1):tt(2);

% Get all good trials
tr_toplot = setdiff(indtrial(D,condition),badtrials(D));

%% Set time points > RT = nan, skip for now

% Get reaction times for subject

trial_data = zeros(length(tr_toplot),2);


for t = 1:length(tr_toplot)

    
    
    trial_data(t,1) = tr_toplot(t);  
    trial_data(t,2) = behav.D.trials(t).events(1).duration;  
    
    
end

% Sort Matrix
[trial_data(:,2),idx] = sort(trial_data(:,2));
trial_data(:,1) = trial_data(idx,1);

% Create matrix (trials vs. time points)
data = squeeze(D(elec,idt,trial_data(:,1)));
    
for t = 1:length(tr_toplot)

    inx = D.time(idt) >  trial_data(t,2);
    data(inx,t) = NaN;
        
    
end


%% Plot heat map with imagesc

x = round(D.time(idt)*1000);
y = 1:length(tr_toplot);

data = data';

z = data; 


% if interp
% 
%     %// Define integer grid of coordinates for the above data
%     [X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
% 
%     %// Define a finer grid of points
%     [X2,Y2] = meshgrid(1:0.1:size(data,2), 1:0.1:size(data,1));
% 
%     %// Interpolate the data and show the output
%     outData = interp2(X, Y, data, X2, Y2, 'linear');
%     figure;
%     imagesc(outData);
%     
% end


figure;
%subplot(2,1,1);
pcolor(x,y,z)

if interp
    
    shading 'interp'

else 
    
    shading 'flat'

end



set(gca,'Ydir','reverse')

%subplot(2,1,2);
%pcolor([data nan(nr,1); nan(1,nc+1)]);
% shading flat;

% Set color map
%colormap(flipud(brewermap([],'RdYlBu')))
colormap(cm_new)


h = colorbar

% Set color map boundries
maxb = median(max(z)); 
minb = median(min(z));

% 
% if deactive
%     
%     caxis([minb abs(minb)])
% 
% else
%     
%     caxis([-maxb maxb])
%     
% end


if deactive
    
    caxis([-2 2])

else
    
    caxis([-2 2])
    
end

%% Plotting Properties

hold on
plot([0 0],ylim,'-k','linewidth',2)

%ticker properties
xlabel('Time (ms)','fontsize',16);
ylabel('Math Trials','fontsize',16);
title('SPL: Math-Active','fontsize',16);

 set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'LineWidth'   , 1         );
   set(gca,'linewidth',2)


set(gca,'XTick',min(x):200:max(x))

print('-opengl','-r300','-dpng','Mathdeactive');

