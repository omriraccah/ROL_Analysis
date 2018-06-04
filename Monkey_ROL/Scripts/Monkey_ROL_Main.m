function ROL = Monkey_ROL_Main(fname,channels,condition, bs_frac, deactive_arr)

% channels: labels of channels to look at. (ALL TOGETHER!)
% condition: what category to look at
% bs_frac: what fraction of the data you would like to use for each
% iteration of bootstrapping (default: 1, or %100)
% deactive_arr: an array that hold 1s and 0s for each relevant channel - 1
% = deactivation, 0 = activation
%% Bootrapping Sequence

% Suppress figure output
set(0,'DefaultFigureVisible','on')

% Create folder to save ROL results
if ~exist([pwd, 'ROLs'],'dir')
            mkdir(pwd,'ROLs');
end

%  If empty, select subject file
if isempty(fname)
    fname = spm_select;
end
D = spm_eeg_load(fname);

if isempty(channels)
    channels = chanlabels(D);
end

indchan = indchannel(D,channels);

% Default bootstrapping parameters
nperm = 1000;
nboots = 1000;

if isempty(condition)
    condition = D.condlist;
end

% Get trial condition indexes 
indtr =  setdiff(indtrial(D,condition),badtrials(D));

% Initialize matrix to collect bootstrap estimates
onsets = zeros(nboots,length(indchan));

 % Get onsets for each bootstrap
boot = zeros(ceil(bs_frac*length(indtr)),1);
for ib = 1:nboots
    for it = 1:ceil(bs_frac*length(indtr))
        tmp =randperm(length(indtr));
        boot(it) = tmp(1);
    end
    
    for ic = 1:length(indchan)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                
        % Run ROL Analysis
        Onsdata = Fit_ROL_WideBins(D,indchan(ic),indtr(boot),[],deactive_arr(ic),0);
        onsets(ib,ic) = Onsdata.all;
        
    end
end

%% Perform permutation for each pair of channels

nel = 1;

if nel > 1

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
%     % Save results
    %path = spm_fileparts(D.fname);
     ROL = struct;
     ROL.onsets = onsets;
    ROL.pairs = pair;
    ROL.pOnset_Pair = pOnset_Pair;
    ROL.crit_p = crit_p;
    ROL.hPair = hPair;

    mkdir(strcat(pwd,'/ROLs/', strcat('ROL_results_',fname(end-11: end-4))));

    save(fullfile(pwd,'/ROLs/',strcat('ROL_results_',fname(end-11: end-4)),'ROL'))

else
    
    ROL = struct;
    ROL.onsets = onsets;

    mkdir(strcat(pwd,'/ROLs/', strcat('ROL_results_',fname(end-11: end-4))));

   save(fullfile(pwd,'/ROLs/',strcat('ROL_results_',fname(end-11: end-4)),'ROL'))



end

%% Plot common elec signals

num_active = sum(deactive_arr(:)==0);
num_deactive = sum(deactive_arr(:)==1);


% Get plot colors
active_colors = [0 0 1];
active_colors = repmat(active_colors,num_active,1)


deactive_colors = [8 180 238] ./ 255;
deactive_colors = repmat(deactive_colors,num_deactive,1)

%Plot_Common_Elec_Signals(fname,indchan,{condition},ROL,[],[],[]); % for paper: vertcat(active_colors,deactive_colors)% [1 0 0;.6 0 0] if you want red and light red
%Plot_Common_Elec_Normalized_Signals(fname,indchan,{condition},ROL,[],[],[],[],deactive_arr);


%% Plot singal electrodes
for cp = 1:length(indchan)
    
    %LBCN_plot_averaged_signal_epochs(fname,indchan(cp),{condition},[],[],[],[],[1 0 0]);
    Plot_ROL_Single_Channel(fname,indchan(cp),{condition},ROL,cp,[]);
    
end

%% Get ROL of Average signal & Plot Bar Graph for each channel, respectively
% 
% % Initialize matrix to collect bootstrap estimates
% mean_onsets = zeros(1,length(indchan));
% bin_data = cell(1,length(indchan));
% ROL_Bin_Index = cell(1,length(indchan));
% sig_thresh = cell(1,length(indchan));
% 
% % This will run the ROL for average signals 
% for ic = 1:length(indchan)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
%                 
%         % Run ROL Analysis
%         mean_Onsdata = Fit_ROL(D,indchan(ic),indtr,[],deactive_arr(ic),1);
%         mean_onsets(1,ic) = mean_Onsdata.all;
%         bin_data(1,ic) = {mean_Onsdata.bin_data};
%         ROL_Bin_Index(1,ic) = mean_Onsdata.ROL_Bin_Index;
%         
%         if deactive_arr(ic) == 1
%             
%              sig_thresh(1,ic) = {mean_Onsdata.neg_thr};
%              
%         else
%             
%             sig_thresh(1,ic) = {mean_Onsdata.pos_thr};
%             
%             
%         end
%         
% end
% 
% 
% % Run ROL_Bar function
% ROL_Bar_Plot(fname,indchan,{condition},mean_onsets,bin_data,ROL_Bin_Index,sig_thresh,[],[],[]) % [1 0 0;.6 0 0] if you want red and light red
% 

