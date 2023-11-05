function ai=avg_infoscore(all_infoscore,trials_selection)

ai=[];
for i=1:size(all_infoscore,1)
    ti=all_infoscore(i,trials_selection(i,:));
    ti=nanmean(ti,2);
    ai(i,1)=nanmean(ti);
end

