% script for running paired permutations tests (in the context of ROL: find
% sections in two arrays that have no NaN values and compare those.
function [res] = script_permutations_paired_2conds(values1,values2,nperm,suffix)

if nargin<3 || isempty(nperm)
    nperm = 1000;
end

if nargin<4 || isempty(suffix)
    suffix = '';
end

if numel(values1)~=numel(values2)
    beep
    disp('For paired tests the vectors should have the same dimensions')
end

% permutation = zeros(numel(values1),nperm);
% Compute permutations 


% truediff = median(values1) - median(values2);
truediff = abs(median(values1 - values2));
pchan_WSR = signrank(values1,values2);
perm = zeros(numel(values1),nperm);
permdiff = zeros(nperm,1);
differences = values1 - values2;
for p = 1:nperm
    indperm = rand(numel(values1),1);
    perm(:,p) =  indperm;
    diff = values1 - values2;
    multp = ones(numel(values1),1);
    multp(perm(:,p)>0.5) = -1;
    diffp = diff .* multp;
    permdiff(p) = abs(median(diffp));
end
permdiff = [permdiff;truediff];
pchan_Cond = length(find(permdiff>=truediff))/(nperm + 1);

% Save results
res.pchan_Cond = pchan_Cond;
res.pchan_WSR = pchan_WSR;
res.values = differences;
save(['Results_Compare_',suffix,'.mat'],'res')