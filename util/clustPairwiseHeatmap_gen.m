function CM_out=clustPairwiseHeatmap_gen(CM1,CM2,perm,flipSign)

CM1(logical(eye(size(CM1)))) = 1;
if flipSign==0
    CM_out=flip(CM1(perm,perm));
else
    CM_out=CM1(perm,perm);
end
