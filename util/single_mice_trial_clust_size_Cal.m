region_sz_6_cir_k={};
region_sz_6_cir_shuf_k={};

load([foldername_multiGeo{6},'\','neuronIndividuals_new.mat']);

for i=2:10
    [~,~,region_sz_6_cir_k{i},~,~,~,region_sz_6_cir_shuf_k{i}]=DBSCAN_region_quantify_022422(group_ori_multiGeo_k{6,2}{i},neuronIndividuals_new,[]);      
end