function [avg_ratemaps,avg_countTimes,avg_cpt]=avg_lap_ratemap_countTime(n_dir1,b_dir1)

%% for each lap, calculate ratemap and countTime
for k=1:length(n_dir1)
    [all_ratemaps{k,1},~,~,countTime{k,1},~,~,~,cpt{k,1}] = calculatingCellSpatialForSingleData_Suoqin(n_dir1{k},b_dir1{k}.position,b_dir1{k}.time,b_dir1{k}.ROI,15,1:size(n_dir1{k}.C,1),0.1*max(n_dir1{k}.C,[],2),'S',[],[],[0 inf],5);
end

%% avg ratemap/countTime across laps together
avg_ratemaps=cell(length(all_ratemaps{1,1}),1);
avg_countTimes=[];
avg_cpt=[];

for k=1:length(n_dir1)
    currCell=all_ratemaps{k,1};
    currCountTime=countTime{k,1};
    currCountTime(isnan(currCountTime))=0;
    
    for k1=1:length(currCell)
        if ~isempty(currCell{k1})
            currCell{k1}(isnan(currCell{k1}))=0;
            if isempty(avg_ratemaps{k1})
                avg_ratemaps{k1}=currCell{k1};
            else
                avg_ratemaps{k1}=avg_ratemaps{k1}+currCell{k1};
            end
        end
    end
    avg_ratemaps{k1}=avg_ratemaps{k1}/length(currCell);

    if isempty(avg_countTimes)
        avg_countTimes=countTime{k,1};
        avg_cpt=cpt{k,1};
    else
        avg_countTimes=avg_countTimes+countTime{k,1};
        avg_cpt=avg_cpt+cpt{k,1};
    end
end
avg_countTimes=avg_countTimes/k;
avg_cpt=avg_cpt/k;