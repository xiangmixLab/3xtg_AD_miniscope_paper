function [avg_region_kmean,avg_region_shuf_kmean,avg_region_ICA_kmean,avg_region_ICA_shuf_kmean]=clustNum_ICA_kmean(foldername,cond)

tic;
for i=1:length(foldername)
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    for j=cond
        for k=2:10
            group_kmean_clustNum{i,j}{k}=cluster_determine_by_suoqin_NMF_firstPeakCoph_022422(neuronIndividuals_new{j},100,10,k);
            [~,~,avg_region_kmean{i,j}{k},~,~,~,avg_region_shuf_kmean{i,j}{k}]=DBSCAN_region_quantify_022422(group_kmean_clustNum{i,j}{k},neuronIndividuals_new,[]);
            
            AssemblyTemplates_clustNum{i,j}{k} = assembly_patterns_num(neuronIndividuals_new{j}.C,max(group_kmean_clustNum{i,j}{k}));
            [time_projection_clustNum{i,j}{k}] = assembly_activity(AssemblyTemplates_clustNum{i,j}{k},neuronIndividuals_new{j}.C);

            group_ICA_clustNum{i,j}{k}=AssemblyTemplateCellGroupInfer(AssemblyTemplates_clustNum{i,j}{k},time_projection_clustNum{i,j}{k},neuronIndividuals_new{j}.C,[]);  
            [~,group_ICA_clustNum{i,j}{k},~]=alignClusterIdx(group_kmean_clustNum{i,j}{k},group_ICA_clustNum{i,j}{k});

            [d1,d2]=size(neuronIndividuals_new{1}.Cn);
            [~,~,avg_region_ICA_kmean{i,j}{k},~,~,~,avg_region_ICA_shuf_kmean{i,j}{k}]=DBSCAN_region_quantify_022422(group_kmean_clustNum{i,j}{k},neuronIndividuals_new,[]);

            toc;
        end
    end
end
