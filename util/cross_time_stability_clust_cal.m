% cross_time_stability_arranged
function [gp_aligned,A_color,A_color_region]=cross_time_stability_clust_cal(foldername,gp_trial,trial_to_do)

% 1. calculate cluster
gp_aligned={};
neuron_all={};
nA={};
for tk=1:length(foldername)
    load([foldername{tk},'\','neuronIndividuals_new.mat']);
    for j=trial_to_do
        max_clust=max(gp_trial{tk,j});
        neuronCurr=neuronIndividuals_new{j};
        neuronIndividuals_cross_time=neuronIndividuals_new_split(neuronIndividuals_new{j},15*60*6,(15*60*6)/2);
        for k=1:length(neuronIndividuals_cross_time)
            [~,gp_aligned{tk,j}{k}]=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_cross_time{k},100,10,max_clust);
            neuron_all{tk,j}{k}=neuronIndividuals_cross_time{k};
            neuron_all{tk,j}{k}.Cn=neuronIndividuals_new{j}.Cn;
        end
        
    end
    nA{tk}=neuronIndividuals_new{1}.A;
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
            neuron_all{tk,j}{k}.A=nA{tk};
            [A_color{tk,j}{k},A_color_region{tk,j}{k}]=DBSCAN_region_quantify_022422(gp_aligned{tk,j}{k},{neuron_all{tk,j}{k}},[]);
        end
    end
end
