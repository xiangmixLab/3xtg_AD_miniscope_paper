function [all_ratemaps_1,all_countTimes_1,all_infoscore_1,trial_type]=ratemap_countTime_infoscore_rearrange(all_ratemaps,all_countTimes,all_infoscore,trial_idx)

all_ratemaps_1={};
all_countTimes_1={};
trial_type=[];

ctt=1;
for i=1:size(all_ratemaps,1)
    for j=1:size(all_ratemaps,2)
        for k=1:length(all_ratemaps{i,j})
            all_ratemaps_1{ctt}=all_ratemaps{i,j}{k};
            all_countTimes_1{ctt}=all_countTimes{i,j};
            all_countTimes_1{ctt}(isnan(all_countTimes_1{ctt}(:)))=0;
            
            for tps=1:size(trial_idx,1)
                if ismember(j,trial_idx(tps,:))
                    trial_type(ctt)=tps;
                    break;
                end
            end

            ctt=ctt+1;
        end
    end
end


%% ratemaps with high/low information score
all_infoscore_1=cell2mat(reshape(all_infoscore,[],1));
all_infoscore_1=all_infoscore_1(:,end);