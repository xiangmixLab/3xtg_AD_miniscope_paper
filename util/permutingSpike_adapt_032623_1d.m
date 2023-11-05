
function [place_cells,all_infoScore,all_infoScore_norm,all_coherence] = permutingSpike_adapt_032623_1d(neuron,behavpos,behavtime,temp,occThresh,binsize,small_velo,infosign,coh_thresh)

if ~exist('occThresh','var') || isempty(occThresh)
%     occThresh = 0.2; %061219 adjust code
      occThresh = 0.1;
end

if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end

if ~exist('binsize','var') || isempty(binsize)
    binsize=15;
end

% trim low fr neurons, they tend to have high infoscore but that doesn't
% necessarily mean they are place cells...

%% prepare statistics
nboot = 1000;
countTimeThresh = [0 inf]; 

maxbehavROI=[0 0 max(behavpos(:,1)),max(behavpos(:,2))];

thresh=3*std(neuron.S,[],2);

%% null score
[firingrateAll,countAll,~,countTime,~,~,bininfo,cpt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh,small_velo);

for i=1:length(firingrateAll)
    firingrateAll{i}=nanmean(firingrateAll{i},1);
end
countTime=nanmean(countTime,1);

infoPerSecondnull = zeros(length(firingrateAll),1);
infoPerSpikenull = infoPerSecondnull;

for j = 1:length(firingrateAll)
    MeanFiringRateAll= sum(sum(countAll{j}))/sum(sum(countTime));
    if isempty(firingrateAll{j})
        continue;
    end    
    frt=firingrateAll{j};
    [infoPerSecondnull(j), infoPerSpikenull(j),pxO] = Doug_spatialInfo_020921(frt,cpt,occThresh,MeanFiringRateAll);  
end

%% null coherence

coherencenull = zeros(length(firingrateAll),1);

for j = 1:length(firingrateAll)
    frt=firingrateAll{j};
    [~,coherencenull(j)] = spatial_coherence(frt);  
end

%% shuffled score
S0 = neuron.S;
C0 = neuron.C;
time0 = neuron.time;

infoPerSecondboot = zeros(length(firingrateAll),nboot);
infoPerSpikeboot = infoPerSecondboot;
coherenceboot = infoPerSecondboot;

rng default % for reproduction
firingrateAllt_collection={};
% deltaT=[0.01:(0.99-0.01)/nboot:0.99];

parfor nE = 1:nboot
    S1=trunk_shuffle_data(S0,size(S0,2)); % heavier shuffling, this induce more neurons
    C1=trunk_shuffle_data(C0,size(C0,2));
    neuronboot=[];
    neuronboot.S = S1;
    neuronboot.C = C1;
    neuronboot.time = time0;
    [firingrateAllt,countAllt] = calculatingCellSpatialForSingleData_040321(neuronboot,behavpos,behavtime,maxbehavROI,binsize,1:size(neuronboot.C,1),thresh,temp,[],[],countTimeThresh,small_velo);
    
    for i=1:length(firingrateAllt)
        firingrateAllt{i}=nanmean(firingrateAllt{i},1);
    end

    infoPerSecondbootT = zeros(length(firingrateAllt),1);
    infoPerSpikebootT = zeros(length(firingrateAllt),1);
    coherencebootT = zeros(length(firingrateAllt),1);

    for j = 1:length(firingrateAllt)
        if isempty(firingrateAll{j})
            continue;
        end   
        MeanFiringRateAll= sum(sum(countAllt{1,j}))/sum(sum(countTime));
        frt=firingrateAllt{j};
        [infoPerSecondbootT(j), infoPerSpikebootT(j)] = Doug_spatialInfo_020921(frt, countTime,occThresh,MeanFiringRateAll);   
%         [~,coherencebootT(j)] = spatial_coherence(frt);  
    end

    infoPerSecondboot(:,nE) = infoPerSecondbootT; 
    infoPerSpikeboot(:,nE) = infoPerSpikebootT;
%     coherenceboot(:,nE) = coherencebootT;
    firingrateAllt_collection{nE}=firingrateAllt;
end

%% output prepare
infoScoreSecondboot = [infoPerSecondnull,infoPerSecondboot];

infoScoreSpikeboot = [infoPerSpikenull,infoPerSpikeboot];
coherenceThresh=coh_thresh;
if isequal(infosign,'sec')
    infoScore = infoScoreSecondboot;    
    infoScoreThresh = quantile(infoScore,0.95,2);
    TinfoPerSecond = table([1:length(infoPerSecondnull)]',infoPerSecondnull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);
    
%     coherenceThresh=quantile(coherenceboot,0.95,2);
    place_cells=intersect(place_cells,find(coherencenull>coherenceThresh));
end
if isequal(infosign,'spk')
    infoScore = infoScoreSpikeboot;
    infoScoreThresh = quantile(infoScore,0.95,2);
    TinfoPerSecond = table([1:length(infoPerSpikenull)]',infoPerSpikenull,infoScoreThresh,'VariableNames',{'neuron','infoScore','thresh'});
    TinfoPerSecond2 = sortrows(TinfoPerSecond,{'infoScore'},{'descend'});
    place_cells = TinfoPerSecond2.neuron(TinfoPerSecond2.infoScore > TinfoPerSecond2.thresh);

%     coherenceThresh=quantile(coherenceboot,0.95,2);
    place_cells=intersect(place_cells,find(coherencenull>coherenceThresh));
end
if isequal(infosign,'all')
    infoScore1 = infoScoreSecondboot;      
    infoScore2 = infoScoreSpikeboot; 
    infoScoreThresh1 = quantile(infoScore1,0.95,2);
    infoScoreThresh2 = quantile(infoScore2,0.95,2);
    place_cells = {find(infoPerSecondnull>infoScoreThresh1),find(infoPerSpikenull>infoScoreThresh2)};

%     coherenceThresh=quantile(coherenceboot,0.95,2);
    place_cells={intersect(place_cells{1},find(coherencenull>coherenceThresh)),intersect(place_cells{2},find(coherencenull>coherenceThresh))};

end


all_infoScore=[infoPerSecondnull,infoPerSpikenull];
all_coherence=[coherencenull];
all_infoScore_norm=[(infoPerSecondnull-nanmean(infoScoreSecondboot,2))./nanstd(infoScoreSecondboot,[],2),(infoPerSpikenull-nanmean(infoScoreSpikeboot,2))./nanstd(infoScoreSpikeboot,[],2)];