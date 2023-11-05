function [avg_corr_mice]=corr_border_processing_perMice(all_corr)

avg_corr_mice=[];
for i=1:size(all_corr,1)
    for j=1:size(all_corr,2)
        avg_corr_mice(i,j)=mean(all_corr{i,j});
    end
end
avg_corr_mice=mean(avg_corr_mice,2);




