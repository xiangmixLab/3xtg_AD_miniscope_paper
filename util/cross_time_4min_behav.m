% cross_time_stability_arranged
function [behav_4min]=cross_time_4min_behav(foldername,trial_to_do)

% 1. calculate cluster
behav_4min={};
for tk=1:length(foldername)
    load([foldername{tk},'\','behav.mat']);
    for j=trial_to_do
        max_clust=max(gp_trial{tk,j});
        neuronCurr=neuronIndividuals_new{j};
        neuronIndividuals_cross_time=neuronIndividuals_new_split(neuronIndividuals_new{j},15*60*4);
        for k=1:length(neuronIndividuals_cross_time)
            [~,gp_aligned{tk,j}{k}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_cross_time{k},100,10,max_clust);
            neuron_all{tk,j}{k}=neuronIndividuals_cross_time{k};
        end
    end
end

