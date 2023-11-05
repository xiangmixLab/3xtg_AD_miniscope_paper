function [nnew_half_temp,behavpos_new_half_temp,behavtime_new_half_temp]=OF_cut_half_temporal(neuron,behavpos,behavtime)

%% get the maximal position of x and y axis, to determine how the track is oriented
max_behavpos_x=max(behavpos(:,1));
max_behavpos_y=max(behavpos(:,2));

nnew=neuron;
ntime = double(nnew.time);

behavpos_new_half_temp{1}=behavpos(1:round(size(behavpos,1)/2),:);
behavpos_new_half_temp{2}=behavpos(round(size(behavpos,1)/2):end,:);
behavtime_new_half_temp{1}=behavtime(1:round(size(behavpos,1)/2),:);
behavtime_new_half_temp{2}=behavtime(round(size(behavpos,1)/2):end,:);

nnew_half_temp{1}.C=nnew.C(:,1:round(size(nnew.C,2)/2));
nnew_half_temp{1}.S=nnew.S(:,1:round(size(nnew.C,2)/2));
nnew_half_temp{1}.time=nnew.time(1:round(length(nnew.time)/2));

nnew_half_temp{2}.C=nnew.C(:,round(size(nnew.C,2)/2):end);
nnew_half_temp{2}.S=nnew.S(:,round(size(nnew.C,2)/2):end);
nnew_half_temp{2}.time=nnew.time(round(length(nnew.time)/2):end);

