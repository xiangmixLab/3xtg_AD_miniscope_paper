
for i=1:6
    load([foldername_multiGeo{i},'\','neuronIndividuals_new.mat']);
    avg_region_multiGeo_k{i,3}={};
    avg_region_shuf_multiGeo_k{i,3}={};
    for k=2:10
        [~,~,avg_region_multiGeo_k{i,3}{k},~,~,~,avg_region_shuf_multiGeo_k{i,3}{k}]=DBSCAN_region_quantify_022422(group_ori_multiGeo_k{i,3}{k},neuronIndividuals_new,[]);  
    end
end
