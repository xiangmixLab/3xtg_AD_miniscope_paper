%% manual set parameters
load('F:\AD_square_circle_results_092320\exp_info_etgAD_031621.mat');
AD_idx=[4 4 4 4 3 3 3 3 4 4 3 3 4 4 2 2 2 1 1 1 2 2 1 2 1 1 1 2 2 1 1 1];
folderName=unique(destination);

bad_mice=[7 25];
tk_range=[1:length(folderName)];
tk_range(bad_mice)=[];

behavname={
    'circle1_behav.mat';
    'square1_behav.mat';
    'square2_behav.mat';
    'circle2_behav.mat';
    'circle3_behav.mat';
    'square3_behav.mat';
    'square4_behav.mat';
    'circle4_behav.mat';
    } % the order of behav files may not follow this one

cond_labels={
    'Circle';
    'Square';
    'Square';
    'Circle';
    'Circle';
    'Square';
    'Square';
    'Circle';    
    };

%% load all neurons, behav
all_neuron={};
all_neuron_simp={};
all_keep_idx={};
tic;
for i=1:length(folderName)
    load([folderName{i},'\','neuronIndividuals_new.mat']);
    for j=1:length(behavname)
        all_neuron{i,j}=neuronIndividuals_new{j}.copy;
        
        all_neuron_simp{i,j}.A=all_neuron{i,j}.A;
        all_neuron_simp{i,j}.C=all_neuron{i,j}.C;
        all_neuron_simp{i,j}.C_raw=all_neuron{i,j}.C_raw;
        all_neuron_simp{i,j}.S=C_to_peakS(all_neuron{i,j}.C);
        all_neuron_simp{i,j}.time=all_neuron{i,j}.time;
        all_neuron_simp{i,j}.Cn=all_neuron{i,j}.Cn;
        all_neuron_simp{i,j}.imageSize=all_neuron{i,j}.imageSize;
        all_neuron_simp{i,j}.centroid=all_neuron{i,j}.centroid;
        all_neuron{i,j}=[];
    end
end
toc;

clear all_neuron
% all_behav
all_behav={};
tic;
for i=1:length(folderName)
    for j=1:length(behavname)
        load([folderName{i},'\',behavname{j}]);
        all_behav{i,j}=behav;
        all_behav{i,j}.VidObj=[];
    end
end
toc;

save('all_neuron_behav_050321.mat','all_neuron_simp','all_behav','-v7.3');

% all neuron num
for i=1:size(all_neuron_simp,1)
    all_neuron_num(i,1)=size(all_neuron_simp{i,1}.C,1);
end

%% calculate responsive cells
load('F:\AD_square_circle_results_092320\velocity_041321.mat');
all_neuron_run={};
all_neuron_nrun={};

velo_thresh=0;
thresh={};
for tk=1:32
    for j=1:size(all_neuron_simp,2)
        all_neuron_run{tk,j}=all_neuron_simp{tk,j};
        all_velo{tk,j}=all_velo{tk,j}(1:size(all_neuron_simp{tk,j}.C,2));
        all_neuron_run{tk,j}.C=all_neuron_run{tk,j}.C(:,all_velo{tk,j}>velo_thresh);
        all_neuron_nrun{tk,j}=all_neuron_simp{tk,j};
        all_neuron_nrun{tk,j}.C=all_neuron_nrun{tk,j}.C(:,all_velo{tk,j}<=velo_thresh);
        thresh{tk,j}=0.1*max(all_neuron_simp{tk,j}.C,[],2);
    end
end

all_fr={};
ctt=1;
for tk=1:30    
    [~,~,all_fr{tk,1}]=general_fr_amp_calculation(all_neuron_simp(tk,:),[],[],cell2mat(thresh(tk,:))); % neuron trim done here 0.02Hz can be applied locally
end
all_fr_m=cell2mat(all_fr);
responsive_thresh=0;
% responsive_thresh=quantile(all_fr_m(:),0.25)
responsive_idx={};
for tk=1:30
    for j=1:size(all_fr{tk,1},2)
        responsive_idx{tk,j}=all_fr{tk,1}(:,j)>responsive_thresh;
    end
end

%% modify del_idx, add all neurons that may have infoscore smaller than 0.02
% in any trials
load('F:\AD_square_circle_results_092320\manual_temporal_check_032721.mat')
del_idx=del_ind;

del_idx_1={};
responsive_thresh=0;
for tk=1:32
    [~,~,fr_chk]=general_fr_amp_calculation(all_neuron_simp(tk,:),[],[]); % neuron trim done here
    for i=1:size(fr_chk,2)
        del_idx_1{tk,i}=unique([del_idx{tk}';find(fr_chk(:,i)<=responsive_thresh)]);
    end
end

%% trials selection
trials_selection=[
    1 2 5 6
    2 5 6 8
    2 5 6 8
    1 2 6 8
    2 4 5 6
    2 5 6 8
    5 6 7 8
    1 2 3 5
    3 4 6 8
    2 4 6 8
    1 3 6 8
    1 2 4 7
    2 4 5 6
    2 3 4 5
    1 2 4 6 %% 15
    2 5 6 8
    1 2 4 6 
    1 2 6 8
    1 3 6 8
    1 2 4 6
    3 4 6 8
    1 3 4 6
    2 5 6 8
    2 4 6 8
    1 2 3 4
    1 2 5 6
    1 3 5 7
    1 2 7 8
    1 2 3 4
    1 2 3 5
    1 2 3 4
    1 3 4 7
    ]
    
%% calculate spatial corherence
tic;
binsize=10;
all_spatial_fr={};
all_spatial_count={};
all_spatial_ct={};
all_spatial_coherence={};
for tk=tk_range
%     load([folderName{tk},'\','neuronIndividuals_new.mat']);
    for j=1:size(all_neuron_simp,2)
        [all_spatial_fr{tk,j},all_spatial_count{tk,j},~,all_spatial_ct{tk,j}] = calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{tk,j},all_behav{tk,j}.position,all_behav{tk,j}.time,[0 0 max(all_behav{tk,j}.position(:,1)),max(all_behav{tk,j}.position(:,2))],binsize,1:size(all_neuron_simp{tk,j}.C,1),0.1*max(all_neuron_simp{tk,j}.C,[],2),'S',[],[],[0.1 inf],5);
        for k=1:length(all_spatial_fr{tk,j})
            all_spatial_coherence{tk,j}(k,1)=spatial_coherence(all_spatial_fr{tk,j}{k});
        end
    end
    toc;
end


%% calculate infoscore. no need to use all_neuron_run as the thres is
% applied in the infoscore cal function
% load('all_neuron_behav_041621.mat');
all_pc=cell(32,8);
all_infoscore=cell(32,8);
all_infoscore_norm=cell(32,8);
all_coherence=cell(32,8);
tic;
for i=1:32*8
    try
        [place_cells,infoScore,infoScore_norm,coherencee] = permutingSpike_adapt_040821(all_neuron_simp{i},all_behav{i}.position,all_behav{i}.time,'S',0,10,5,'all',0.3);  

        all_pc{i}=place_cells;
        all_infoscore{i}=infoScore;
        all_infoscore_norm{i}=infoScore_norm;
        all_coherence{i}=coherencee;
    catch
    end
    toc;
end

% equal firing rate sample place cell
save('all_pc_infoscore_050321.mat','all_pc','all_infoscore','all_infoscore_norm','-v7.3');

%% field size
all_fr_spatial=cell(30,8);
all_field_size=cell(30,8);

tic;
for i=1:240
    try
        field_size=[];
        
        [firingrateAllt,countAllt,~,countTimeAllt] = calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{i},all_behav{i}.position,all_behav{i}.time,all_behav{i}.ROI,10,1:size(all_neuron_simp{i}.C,1),0.1*max(all_neuron_simp{i}.C,[],2),'S',[],[],[0.1 inf],5);
        all_fr_spatial{i}=firingrateAllt;
        for j=1:length(firingrateAllt)
            
            if ~isempty(firingrateAllt{j}) && sum(firingrateAllt{j}(:))~=0
                fr=filter2DMatrices_021521(firingrateAllt{j},5,2);  
                fr=fr{1};
                fr=fr>0.5*max(fr(:));
                stats=regionprops(fr);
                field_size(j)=max([stats.Area]);
            else
                field_size(j)=nan;
            end
        end
        
        all_field_size{i}=field_size';
    catch
    end
end
toc;

%% sparsity
all_sparsity=cell(32,8);

tic;
for i=1:32*8
    try
        spr=[];
        
        [firingrateAllt,countAllt,~,countTimeAllt,~,~,~,cpt] = calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{i},all_behav{i}.position,all_behav{i}.time,all_behav{i}.ROI,10,1:size(all_neuron_simp{i}.C,1),0.1*max(all_neuron_simp{i}.C,[],2),'S',[],[],[0.1 inf],5);
        Po=cpt/sum(cpt(:));
        for j=1:length(firingrateAllt)
            
            if ~isempty(firingrateAllt{j}) && sum(firingrateAllt{j}(:))~=0
                spr(j)=sum(sum(firingrateAllt{j}.*Po))^2/sum(sum(firingrateAllt{j}.^2.*Po));
            else
                spr(j)=nan;
            end
        end
        
        all_sparsity{i}=spr';
    catch
    end
end
toc;

