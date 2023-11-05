function [all_sparsity_dir1,all_sparsity_dir2]=lap_sparsity_cal(n_dir1,n_dir2,b_dir1,b_dir2)

all_sparsity_dir1={};
all_sparsity_dir2={};

% ratemaps
for i=1:size(n_dir1,1)
    for j=1:size(n_dir1,2)

        % dir1
        [avg_ratemaps,~,avg_cpt]=avg_lap_ratemap_countTime(n_dir1{i,j},b_dir1{i,j});

        all_sparsity_dir1_temp=[];
        Po=avg_cpt/sum(avg_cpt(:));
        for k1=1:length(avg_ratemaps)
            if ~isempty(avg_ratemaps{k1})
                all_sparsity_dir1_temp(k1,1)=sum(sum(avg_ratemaps{k1}.*Po))^2/sum(sum(avg_ratemaps{k1}.^2.*Po));
            else
                all_sparsity_dir1_temp(k1,1)=nan;
            end
        end
        all_sparsity_dir1{i,j}=all_sparsity_dir1_temp;

        % dir2
        [avg_ratemaps,~,avg_cpt]=avg_lap_ratemap_countTime(n_dir2{i,j},b_dir2{i,j});

        
        Po=avg_cpt/sum(avg_cpt(:));
        for k1=1:length(avg_ratemaps)
            all_sparsity_dir2_temp=[];
            if ~isempty(avg_ratemaps{k1})
                all_sparsity_dir2_temp(k1,1)=sum(sum(avg_ratemaps{k1}.*Po))^2/sum(sum(avg_ratemaps{k1}.^2.*Po));
            else
                all_sparsity_dir1_temp(k1,1)=nan;
            end
        end
        all_sparsity_dir2{i,j}=all_sparsity_dir2_temp;
    end
end
