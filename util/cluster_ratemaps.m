function [ensembleMap,firingrateAll,countTime]=cluster_ratemaps(neuron,idx,behav,temp)

neuron.C=neuron.C(idx,:);
neuron.S=neuron.S(idx,:);

behavpos=behav.position;
behavtime=behav.time;
maxbehavROI=[0,0,max(behavpos,[],1)];
binsize=10;
countTimeThresh=[0 inf];

[firingrateAll,countAll,~,countTime] = calculatingCellSpatialForSingleData_Suoqin(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),[],temp,[],[],countTimeThresh,0);

%% avg ensemble map
ensembleMap=firingrateAll{1};
mapNum=1;
for i=2:length(firingrateAll)
    if ~isempty(firingrateAll{i})
        ensembleMap=ensembleMap+filter2DMatrices(firingrateAll{i},1);
        mapNum=mapNum+1;
    end
end

ensembleMap=ensembleMap/mapNum;

%% trace ensemble map
% n_ensem=Sources2D;
% n_ensem.C=mean(neuron.C,1);
% n_ensem.S=mean(neuron.S,1);
% n_ensem.time=neuron.time;
% 
% [ensemMap] = calculatingCellSpatialForSingleData_040321(n_ensem,behavpos,behavtime,maxbehavROI,binsize,1,3*std(n_ensem.S,[],2),temp,[],[],countTimeThresh,small_velo);
% ensembleMap=ensemMap{1};

