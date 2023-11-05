%% C: CNMF-E footprints
ctt=1;
mname={'Ntg young','AD young','Ntg old','AD old'};
for i=[5 3 20 21]
    load([folderName{i},'\','all_neurons.mat']);
    neurons{1}.A(:,del_ind{i})=[];
    footprint=spatial_footprint_calculation(neurons{1},0.5);
    subplot(1,4,ctt);
    imagesc(footprint);
    title([mname{ctt},' ',num2str(size(neurons{1}.C,1))]);
    ctt=ctt+1;
end

for i=1:32
    load([folderName{i},'\','all_neurons.mat']);
    num_neuron=[];
    for j=1:length(neurons)
        num_neuron(j)=size(neurons{j}.C,1);
    end
    
    all_neuron_num_perSession_raw{i,1}=num_neuron;
    all_neuron_num_perSession(i,1)=round(mean(num_neuron));
end

cell_num=fill_nan_to_cellmat({all_neuron_num_perSession(AD_idx==3),all_neuron_num_perSession(AD_idx==4),all_neuron_num_perSession(AD_idx==1),all_neuron_num_perSession(AD_idx==2)})

%% ratemap examples
bad_mice=[7];
tk_range=[1:length(folderName)];
tk_range(bad_mice)=[];

load('F:\AD_square_circle_results_092320\all_neuron_behav_032921.mat');
for i=tk_range
    for j=1:8
        [fr{i,j},ct{i,j},~,ctime{i,j}]=calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{i,j},all_behav{i,j}.position,all_behav{i,j}.time,all_behav{i,j}.ROI,10,1:size(all_neuron_simp{i,j}.C,1),0.1*max(all_neuron_simp{i,j}.C,[],2),'S',[],[],[0.1 inf],5);
    end
end

% fr amp select


% circle
figure;
ctt=1;
ctt1=1;


midx=[5 13 30 16];
trial_sle=[4 1 1 5];
cell_idx=[62 66 5 15];

for i=midx
subplot(2,6,ctt)
plot(all_neuron_simp{i,trial_sle(ctt1)}.C(cell_idx(ctt1),:));
hold on;
plot(1:size(all_neuron_simp{i,trial_sle(ctt1)}.C,2),ones(1,size(all_neuron_simp{i,trial_sle(ctt1)}.C,2))*0.1*max(all_neuron_simp{i,trial_sle(ctt1)}.C(cell_idx(ctt1),:)),'--','color','k','lineWidth',2);
ylim([0 80]);
ctt=ctt+1;
subplot(2,6,ctt)
trajectory_firingpos_plot(all_behav{i,trial_sle(ctt1)},all_neuron_simp{i,trial_sle(ctt1)},all_behav{i,trial_sle(ctt1)}.ROI,cell_idx(ctt1),0.1*max(all_neuron_simp{i,trial_sle(ctt1)}.C(cell_idx(ctt1),:)),10,'S');
ctt=ctt+1;
subplot(2,6,ctt)
ratemap_plot(fr{i,trial_sle(ctt1)}{cell_idx(ctt1)},ctime{i,trial_sle(ctt1)},1,0,[0 0.5]);
ctt=ctt+1;

ctt1=ctt1+1;
end
set(gcf,'renderer','painters');

% square
figure;
ctt=1;
ctt1=1;
trial_sle=[3 2 2 6];
for i=midx
subplot(2,6,ctt)
plot(all_neuron_simp{i,trial_sle(ctt1)}.C(cell_idx(ctt1),:));
hold on;
plot(1:size(all_neuron_simp{i,trial_sle(ctt1)}.C,2),ones(1,size(all_neuron_simp{i,trial_sle(ctt1)}.C,2))*0.1*max(all_neuron_simp{i,trial_sle(ctt1)}.C(cell_idx(ctt1),:)),'--','color','k','lineWidth',2);
ylim([0 80]);
ctt=ctt+1;
subplot(2,6,ctt)
trajectory_firingpos_plot(all_behav{i,trial_sle(ctt1)},all_neuron_simp{i,trial_sle(ctt1)},all_behav{i,trial_sle(ctt1)}.ROI,cell_idx(ctt1),0.1*max(all_neuron_simp{i,trial_sle(ctt1)}.C(cell_idx(ctt1),:)),10,'S');
ctt=ctt+1;
subplot(2,6,ctt)
ratemap_plot(fr{i,trial_sle(ctt1)}{cell_idx(ctt1)},ctime{i,trial_sle(ctt1)},1,0,[0 0.5]);
ctt=ctt+1;

ctt1=ctt1+1;
end
set(gcf,'renderer','painters');
