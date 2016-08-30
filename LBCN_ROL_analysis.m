function [ROL] = LBCN_ROL_analysis(fname,channels,condition)

% fname: name of SPM MEEG object, .mat
% channels: labels of channels to look at. (ALL TOGETHER!)
% condition: what category to look at

if isempty(fname)
    fname = spm_select;
end
D = spm_eeg_load(fname);

if isempty(channels)
    channels = chanlabels(D);
end

indchan = indchannel(D,channels);

%Defaults
nperm = 1000;
nboots = 1000;

if isempty(condition)
    condition = D.condlist;
end

indtr = indtrial(D,condition);

onsets = zeros(nboots,length(indchan));
 % Get onsets for each bootstrap
 boot = zeros(length(indtr),1);
for ib = 1:nboots
    for it = 1:length(indtr)
        tmp =randperm(length(indtr));
        boot(it) = tmp(1);
    end
    
    for ic = 1:length(indchan)
        % Call ROL code to get the onset for that channel
        Onsdata = respfunc_epoch_bootstrap(D,indchan(ic),indtr(boot));
        onsets(ib,ic) = Onsdata.all;
    end
end


%% Perform permutation for each pair of channels
nel = length(indchan);
npairs = factorial(nel) / (factorial(nel-2)*2);
pair = zeros(npairs,2);
cnt = 1;
for i = 1:nel
    for j = (i+1):nel
        pair(cnt,1) = i;
        pair(cnt,2) = j;
        cnt = cnt+1;
    end
end
    
permutation = zeros(nboots,nperm);
pOnset_Pair = zeros(npairs,1);
perm = zeros(nboots,nperm);
permdiff = zeros(nperm,1);
truediff = zeros(npairs,1);
for i = 1:npairs
    % Compute true difference
    onsch1 = onsets(:,pair(i,1));
    onsch2 = onsets(:,pair(i,2));    
    truediff(i) = nanmean(onsch2 - onsch1);
    for p = 1:nperm
        if i == 1
            % Need to set the permutation matrix
            indperm = rand(nboots,1);
            permutation(:,p) = indperm;
        end
        perm(:,p) =  permutation(:,p);
        diff = onsch2 - onsch1;
        multp = ones(nboots,1);
        multp(perm(:,p)>0.5) = -1;
        diffp = diff .* multp;
        permdiff(p) = abs(nanmean(diffp));
    end
    permdiff = [permdiff;abs(truediff(i))];
    pOnset_Pair(i) = length(find(permdiff>=abs(truediff(i))))/(nperm + 1);
end

% Correct for multiple comparisons using FDR
[crit_p,hPair] = LBCN_FDRcorrect(pOnset_Pair);
% Save results
path = spm_fileparts(D.fname);
ROL = struct;
ROL.onsets = onsets;
ROL.pairs = pair;
ROL.pOnset_Pair = pOnset_Pair;
ROL.crit_p = crit_p;
ROL.hPair = hPair;
save(fullfile(path,'ROL_results.mat'),'ROL')

%% Plot results for each channel
for ic = 1:length(indchan)
    data = squeeze(D(indchan(ic),:,indtr));
    plot_onset(onsets(:,ic),data,D,indchan(ic),Onsdata.time,indtr);
end

function [] = plot_onset(all_resp,data,D,elec,time,indcond)


% calulate 95% confidence interval
SE= nanmedian(all_resp)/sqrt(1000); %% Standard Error
ME= (SE) *2 ;       %% margine of Error
CI= ([ME+nanmedian(all_resp) , nanmedian(all_resp)-ME]);

% % signal plot


figure('Position', [80, 80, 600, 600]); set(1,'DefaultFigureVisible','on')
subplot(2,1,1)
plot(time,nanmean(data,2),'k-','LineWidth',3)
hold on
plot([nanmedian(all_resp) nanmedian(all_resp)],ylim,'r-','LineWidth',1.2)

patient_id = D.fname;
fig_name = fullfile(D.path,strcat(patient_id(end-11: end-4),char(chanlabels(D,elec)), '.jpg'));
%line1 = plot([CI(1) CI(1)], ylim,'b')
%line2 = plot([CI(2) CI(2)], ylim,'b')
xlim([time(1) 0.45])
xlabel('Time (s)')
ylabel('Signal')

y = ylim;
title (['t=', num2str(nanmedian(all_resp))]);

% Add a patch
patch([CI(2), CI(1), CI(1), CI(2)],[y(1), y(1), y(2), y(2)],'r','EdgeColor','red','EdgeAlpha','0.1')
alpha(0.1)


% create figure for Zeinab trial-by-trial
subplot(2,1,2)
imagesc(time, 1:size(data,2), data');
hold on
plot ([nanmedian(all_resp) nanmedian(all_resp)], ylim, 'r','LineWidth',2);
xlabel('Time (s)')
xlim([time(1) 0.45])
ylabel('Trials')

% save figure
saveas(gcf, fig_name);
movefile(fig_name, [pwd '/figs']);

% save(strcat('/Users/parvizilab/Downloads/ROL_Scripts/resp/',...
%       'all_resp_', file_name(end-12:end-4), '_',num2str(elec(e))));
%   
% adata= all_resp;
% save(strcat('/Users/parvizilab/Downloads/ROL_Scripts/adata/',...
%       'adata', file_name(end-12:end-4), '_',num2str(elec(e))));
       
    

