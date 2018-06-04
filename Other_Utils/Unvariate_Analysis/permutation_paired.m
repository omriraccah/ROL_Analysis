function p_actual = permutation_paired(adata, bdata, reps)

numa = length(adata);
numb = length(bdata);

real_meandiff = nanmean(adata(:))-nanmean(bdata(:));

permdiff = NaN*ones(1,reps);
for r = 1:reps
    inds = find(rand(1,numa)>0.5);
    tempa = adata;
    tempb = bdata;
    tempa(inds) = bdata(inds);
    tempb(inds) = adata(inds);
    permdiff(r) = nanmean(tempa)-nanmean(tempb);
end

if (real_meandiff>0)
    p = sum(real_meandiff>[real_meandiff permdiff])/(reps+1);
else
    p = sum(real_meandiff>=[real_meandiff permdiff])/(reps+1);
end

p_actual = min(p,1-p)*2;

end
