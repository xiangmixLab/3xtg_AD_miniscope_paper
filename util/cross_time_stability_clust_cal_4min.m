% cross_time_stability_arranged
function [gp_aligned,A_color,A_color_region,behav_all,neuron_all]=cross_time_stability_clust_cal_4min(foldername,gp_trial,trial_to_do)

% 1. calculate cluster
gp_aligned={};
neuron_all={};
behav_all={};
for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    load([foldername{tk},'\','behav.mat']);
    for j=trial_to_do
        max_clust=max(gp_trial{tk,j});
        neuronCurr=neuronIndividuals_new{j};
        neuronIndividuals_cross_time=neuronIndividuals_new_split(neuronIndividuals_new{j},15*60*4,(15*60*4));
        for k=1:length(neuronIndividuals_cross_time)
            [~,gp_aligned{tk,j}{k}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_cross_time{k},100,10,max_clust);
            neuron_all{tk,j}{k}=neuronIndividuals_cross_time{k};
        end
        behav_all{tk,j}=behav_split(behavIndividuals{j},30*60*4,30*60*4);
    end
end

% 2. align clusters
for tk=1:size(gp_aligned,1)
    for j=1:size(gp_aligned,2)
        for k=2:length(gp_aligned{tk,j})
            [~,gp_aligned{tk,j}{k}]=determineSharedCells_new(gp_aligned{tk,j}{k-1},gp_aligned{tk,j}{k});
        end
    end
end

% 3. colored footprints and regions
A_color={};
A_color_region={};
for tk=1:size(gp_aligned,1)
    for j=1:size(gp_aligned,2)
        for k=1:length(gp_aligned{tk,j})
            neuron_all{tk,j}{k}.imageSize=[240,376];
            [A_color{tk,j}{k},A_color_region{tk,j}{k}]=DBSCAN_region_quantify_022422(gp_aligned{tk,j}{k},{neuron_all{tk,j}{k}},[]);
        end
    end
end
