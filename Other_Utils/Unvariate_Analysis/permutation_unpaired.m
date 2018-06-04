function p_actual = permutation_unpaired(adata, bdata, reps)

numa = length(adata);
numb = length(bdata);
numtot = numa+numb;

if (size(adata,1)>size(adata,2))
    data = vertcat(adata, bdata);
else
    data = horzcat(adata, bdata);
end

real_meandiff = nanmean(adata(:))-nanmean(bdata(:));

permdiff = NaN*ones(1,reps);
for r = 1:reps
    inds = randperm(numtot);
    indsa = inds(1:numa);
    indsb = inds(numa+1:end);
    permdiff(r) = nanmean(data(indsa))-nanmean(data(indsb));
end

p = sum(real_meandiff>[permdiff real_meandiff])/(reps+1);

p_actual = min(p,1-p)*2;


end


