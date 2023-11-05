function [fr_all,ctime_all] = calculatinglinearTrackRateMaps_laps_noDir(ndat,bdat,binsize,temp)


for k1=1:size(ndat,1)
    for k2=1:size(ndat,2)
        if ~isempty(ndat{k1,k2})
            for k3=1:size(ndat{k1,k2},1)
                for k4=1:size(ndat{k1,k2},2)
                    if ~isempty(ndat{k1,k2}{k3,k4})
                        [fr_all{k1,k2}{k3,k4},ct_all{k1,k2}{k3,k4},~,ctime_all{k1,k2}{k3,k4}]=calculatingCellSpatialForSingleData_040321(ndat{k1,k2}{k3,k4},bdat{k1,k2}{k3,k4}.position,bdat{k1,k2}{k3,k4}.time,bdat{k1,k2}{k3,k4}.ROI,binsize,1:size(ndat{k1,k2}{k3,k4}.C,1),3*std(ndat{k1,k2}{k3,k4}.S,[],2),temp,[],[],[0 inf],0);
                    end
                end
            end
        end
    end
end
