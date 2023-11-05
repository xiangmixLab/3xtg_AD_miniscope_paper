function [neuron_dir1,neuron_dir2,behav_dir1,behav_dir2]=lap2dir(neuron_lap,behav_lap)

neuron_dir1={};
neuron_dir2={};
behav_dir1={};
behav_dir2={};

for i=1:size(neuron_lap,1)
    for j=1:size(neuron_lap,2)
        neuron_dir1{i,j}=neuron_lap{i,j}(:,1);
        neuron_dir2{i,j}=neuron_lap{i,j}(:,2);
        behav_dir1{i,j}=behav_lap{i,j}(:,1);
        behav_dir2{i,j}=behav_lap{i,j}(:,2);
        
        neuron_dir1{i,j}=neuron_dir1{i,j}(~cellfun('isempty',neuron_dir1{i,j}));
        neuron_dir2{i,j}=neuron_dir2{i,j}(~cellfun('isempty',neuron_dir2{i,j}));
        behav_dir1{i,j}=behav_dir1{i,j}(~cellfun('isempty',behav_dir1{i,j}));
        behav_dir2{i,j}=behav_dir2{i,j}(~cellfun('isempty',behav_dir2{i,j}));
    end
end