function [rateMap_corr_all_halves]=half_trial_ratemap_corr(behavpos_half_ratemap,behavpos_half_ct)



% half trial 
rateMap_corr_all_halves=[];
tic;
for tk=1:size(behavpos_half_ratemap,1)
    rateMap_corr_temp=[];
   for i=1:size(behavpos_half_ratemap{tk,1},2)
       rateMap_corr_temp(:,i)=rateMap_correlation(behavpos_half_ratemap{tk,i}{1},behavpos_half_ratemap{tk,i}{2},behavpos_half_ct{tk,i}{1},behavpos_half_ct{tk,i}{2},1,1);
       % rateMap_corr_temp: cell * trial matrix
   end
   rateMap_corr_all_halves(tk,1)=nanmean(rateMap_corr_temp(:));

toc;
end