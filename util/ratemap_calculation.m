function [half_trial_ratemaps,full_trial_ratemaps,half_trial_countTime,full_trial_countTime]=ratemap_calculation(neuron,behav)

temp='S';

%% full ratemap

behavpos=behav.position;
behavtime=behav.time;
maxbehavROI=[0,0,max(behavpos,[],1)];
binsize=10;
thresh=[];
countTimeThresh=[0.1 inf];

[full_trial_ratemaps,~,~,full_trial_countTime] = calculatingCellSpatialForSingleData_Suoqin(neuron,behavpos,behavtime,maxbehavROI,binsize,1:size(neuron.C,1),thresh,temp,[],[],countTimeThresh,0);

%% half ratemap
half_trial_ratemaps={};
half_trial_countTime={};

% s_n=round(size(neuron.C,2)/2);
% s_b=round(size(behav.position,1)/2);

% first half in time
% n=neuron;
% n.C=neuron.C(:,1:s_n);
% n.S=neuron.S(:,1:s_n);
% n.time=n.time(1:s_n);
% 
% behavpos=behav.position(1:s_b,:);
% behavtime=behav.time(1:s_b,:);

[n,behavposs,behavtimee]=OF_cut_half_temporal(neuron,behavpos,behavtime);

thresh=3*std(n{1}.S,[],2);

[half_trial_ratemaps{1},~,~,half_trial_countTime{1}] = calculatingCellSpatialForSingleData_Suoqin(n{1},behavposs{1},behavtimee{1},maxbehavROI,binsize,1:size(n{1}.C,1),thresh,temp,[],[],countTimeThresh,10);

% second half
% n=neuron;
% n.C=neuron.C(:,s_n+1:end);
% n.S=neuron.S(:,s_n+1:end);
% n.time=n.time(s_n+1:end);
% 
% behavpos=behav.position(s_b+1:end,:);
% behavtime=behav.time(s_b+1:end,:);
thresh=3*std(n{2}.S,[],2);

[half_trial_ratemaps{2},~,~,half_trial_countTime{2}] = calculatingCellSpatialForSingleData_Suoqin(n{2},behavposs{2},behavtimee{2},maxbehavROI,binsize,1:size(n{2}.C,1),thresh,temp,[],[],countTimeThresh,10);
