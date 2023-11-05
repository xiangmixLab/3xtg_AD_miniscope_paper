function [behavpos_half_ratemap,fr_all,behavpos_half_ct,ct_all]=half_and_full_trial_ratemap_calculation(foldername)

behavpos_half_ratemap={};
behavpos_half_ct={};
fr_all={};
ct_all={};

for i=1:length(foldername)
    load([foldername{i},'\','neuronIndividuals_new.mat']);
    load([foldername{i},'\','behav.mat']);
    for j=1:length(behavIndividuals)
        [behavpos_half_ratemap{i,j},fr_all{i,j},behavpos_half_ct{i,j},ct_all{i,j}]=ratemap_calculation(neuronIndividuals_new{j},behavIndividuals{j});
    end
end