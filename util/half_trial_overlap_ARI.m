function [overlap_half,overlap_half_shuf]=half_trial_overlap_ARI(group_half,shuffle)

overlap_half=[];
overlap_half_shuf=[];

for i=1:size(group_half,1)
    for j=1:size(group_half,2)
        [~,~,~,overlap_half(i,j)]=new_cluster_similarity_index_031523(group_half{i,j}{1},group_half{i,j}{2});

        ohft=[];
        for s=1:shuffle
            [~,~,~,ohft(s)]=new_cluster_similarity_index_031523(group_half{i,j}{1}(randperm(length(group_half{i,j}{1}))),group_half{i,j}{2}(randperm(length(group_half{i,j}{2}))));
        end
        overlap_half_shuf(i,j)=nanmean(ohft);
    end
end
overlap_half=nanmean(overlap_half,2);
overlap_half_shuf=nanmean(overlap_half_shuf,2);