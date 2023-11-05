%% figure 4 linear track
load('F:\AD_3xtg_linear_track_013121\dat_info.mat');
binsize=15;
load('F:\AD_3xtg_linear_track_013121\all_neuron_behav_laps_042221.mat')
load('F:\AD_3xtg_linear_track_013121\all_neuron_num_042221.mat')
AD_idx=[4 2 4 2 4 2 4 3 3 1 3 3 4 4 3 1 3 1 4 4 2 2 2 1 1];

folderName=foldername;

bad_mice=[];
tk_range=[1:length(folderName)];
tk_range(bad_mice)=[];

behavname={
    'hor1_behav.mat';
    'hor2_behav.mat';
    'hor3_behav.mat';
    'vet1_behav.mat';
    }

cond_labels={
    'hor1';
    'hor2';
    'hor3';
    'vet1';
    };

trials_selection=repmat([1:4],25,1);

%% panel A-D: ratemap examples
binsize=15;
fr_all={};
ct_all={};
ctime_all={};
for tk=1:size(all_neuron_simp_laps1,1)
    for i=1:size(all_neuron_simp_laps1,2)
        [fr_all{tk,i},ct_all{tk,i},~,ctime_all{tk,i}]=calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp_water_rm{tk,i},all_behav_water_rm{tk,i}.position,all_behav_water_rm{tk,i}.time,all_behav{tk,i}.ROI,binsize,1:size(all_neuron_simp_water_rm{tk,i}.C,1),0.1*max(all_neuron_simp_water_rm{tk,i}.C,[],2),'S',[],[],[0 inf],5);
    end
end

figure;
ctt=1;
for tk=[7 9 18 23]
%     subplot(4,4,ctt);
%     trajectory_firingpos_plot(all_behav_water_rm{tk,1},all_neuron_simp_water_rm{tk,1},all_behav{tk,1}.ROI,10,max(all_neuron_simp_water_rm{tk,1}.C(10,:))*0.1,binsize,'S');
%     
%     subplot(4,4,ctt+4);
%     ratemap_plot(nanmean(fr_all{tk,1}{10},1),nanmedian(ctime_all{tk,1},1),1,1,[]);
%     caxis([0 0.15]);
    subplot(4,2,ctt);
    trajectory_firingpos_plot(all_behav_water_rm{tk,1},all_neuron_simp_water_rm{tk,1},all_behav{tk,1}.ROI,15,max(all_neuron_simp_water_rm{tk,1}.C(15,:))*0.1,binsize,'S');
    subplot(4,2,ctt+4);
    ratemap_plot(nanmean(fr_all{tk,1}{15},1),nanmedian(ctime_all{tk,1},1),1,1,[]);  
    caxis([0 0.15]);
    ctt=ctt+1;
end

figure;
for tk=1:25
    subplot(5,5,tk)
    trajectory_firingpos_plot(all_behav{tk,1},all_neuron_simp{tk,1},all_behav{tk,1}.ROI,10,max(all_neuron_simp{tk,i}.C(10,:))*0.1,binsize,'S');
end

%% E: all neurons, field pos along track, sorted
binsize=15


% ratemap laps
fr_all={};
ct_all={};
ctime_all={};
for tk=1:size(all_neuron_simp_laps1,1)
    for i=1:size(all_neuron_simp_laps1,2)
        ndat=all_neuron_simp_laps1{tk,i};
        bdat=all_behav_laps1{tk,i};
        for k1=1:size(ndat,1)
            for k2=1:size(ndat,2)
                if ~isempty(ndat{k1,k2})
                    [fr_all{tk,i}{k1,k2},ct_all{tk,i}{k1,k2},~,ctime_all{tk,i}{k1,k2}]=calculatingCellSpatialForSingleData_040321(ndat{k1,k2},bdat{k1,k2}.position,bdat{k1,k2}.time,bdat{k1,k2}.ROI,binsize,1:size(ndat{k1,k2}.C,1),0.1*max(ndat{k1,k2}.C,[],2),'S',[],[],[0 inf],0);
                end
            end
        end
    end
end

% ratemap laps avg
[fr_all_dir,ct_all_dir,ctime_all_dir]=laps_avg_fr_ct(fr_all,ct_all,ctime_all,all_neuron_num);

% smooth
fr_all_sm=fr_all_dir;
for tk=1:size(all_neuron_simp_laps1,1)
    for i=1:size(all_neuron_simp_laps1,2)
        for k=1:length(fr_all_sm{1}{tk,i})
            fr_all_sm{1}{tk,i}{k}=filter2DMatrices(fr_all_sm{1}{tk,i}{k},1);
            fr_all_sm{2}{tk,i}{k}=filter2DMatrices(fr_all_sm{2}{tk,i}{k},1);
        end
    end
end

rm_line{1}=linearTrack_ratemap2line(fr_all_sm{1});
rm_line{2}=linearTrack_ratemap2line(fr_all_sm{2});

% place cell select
% for i=1:size(rm_line,1)
%     for j=1:size(rm_line,2)
%         rm_line{i,j}=rm_line{i,j}(all_pc{i,j}{2},:);
%     end
% end
% binarize, merge mice by types 

r_all_laps_avg_bin={};
for tk=1:size(all_neuron_simp_laps1,1)
    for i=1:size(all_neuron_simp_laps1,2)
        fr_all_laps_avg_bin{1}{tk,i}=rm_line{1}{tk,i}>0.99*max(rm_line{1}{tk,i},[],2);
        fr_all_laps_avg_bin{2}{tk,i}=rm_line{2}{tk,i}>0.99*max(rm_line{2}{tk,i},[],2);
    end
end

for i=1:2
    [fr_sm_sorted,fr_sm_ori,fr_sm_hc,idx_listo]=figure6_laps_ratemap_show(fr_all_laps_avg_bin{i},rm_line{i},AD_idx,[])
    % [fr_bin_111]=linemap_apply_sorting(fr_bin_11,idxx1);

    % illustration
    figure;
    for j=[1:4]
        subplot(4,8,2*j-1)
        imagesc(fr_sm_hc{1,j});caxis([0 1.0])
        subplot(4,8,2*j-1+8)
        imagesc(fr_sm_hc{2,j});caxis([0 1.0])
        subplot(4,8,2*j-1+16)
        imagesc(fr_sm_hc{3,j});caxis([0 1.0])
        subplot(4,8,2*j-1+24)
        imagesc(fr_sm_hc{4,j});caxis([0 1.0])
    end
    colormap(jet);
end

%% new E, F: fr, amp, both dir merge
% load('all_neuron_behav_waterrm_041621.mat');
fr_dir1_run={};
fr_dir2_run={};

amp_dir1_run={};
amp_dir2_run={};

del_idx=cell(25,1);
[neuron_dir]=laps_avg_neuroDat(all_neuron_simp_laps1,all_behav_laps1);

for tk=tk_range
%     trial_select_cir=[];
%     trial_select_sqr=[];
%     for j=1:length(trials_selection(tk,:))
%         if ~isempty(strfind(cond_labels{trials_selection(tk,j)},'hor'))
%             trial_select_cir=[trial_select_cir,trials_selection(tk,j)];
%         else
%             trial_select_sqr=[trial_select_sqr,trials_selection(tk,j)];
%         end
%     end
    
    [fr_dir1_run{tk},amp_dir1_run{tk}]=general_fr_amp_calculation(neuron_dir{1}(tk,:),del_idx{tk},0); % neuron trim done here 0.02Hz can be applied locally
    [fr_dir2_run{tk},amp_dir2_run{tk}]=general_fr_amp_calculation(neuron_dir{2}(tk,:),del_idx{tk},0);

end
fr_amp_dat={};
fr_amp_dat{1,1}=cell2mat(fr_dir1_run(:,AD_idx==3)');
fr_amp_dat{1,2}=cell2mat(fr_dir1_run(:,AD_idx==4)');
fr_amp_dat{1,3}=cell2mat(fr_dir1_run(:,AD_idx==1)');
fr_amp_dat{1,4}=cell2mat(fr_dir1_run(:,AD_idx==2)');

fr_amp_dat{2,1}=cell2mat(fr_dir2_run(:,AD_idx==3)');
fr_amp_dat{2,2}=cell2mat(fr_dir2_run(:,AD_idx==4)');
fr_amp_dat{2,3}=cell2mat(fr_dir2_run(:,AD_idx==1)');
fr_amp_dat{2,4}=cell2mat(fr_dir2_run(:,AD_idx==2)');


fr_amp_dat{3,1}=cell2mat(amp_dir1_run(:,AD_idx==3)');
fr_amp_dat{3,2}=cell2mat(amp_dir1_run(:,AD_idx==4)');
fr_amp_dat{3,3}=cell2mat(amp_dir1_run(:,AD_idx==1)');
fr_amp_dat{3,4}=cell2mat(amp_dir1_run(:,AD_idx==2)');


fr_amp_dat{4,1}=cell2mat(amp_dir2_run(:,AD_idx==3)');
fr_amp_dat{4,2}=cell2mat(amp_dir2_run(:,AD_idx==4)');
fr_amp_dat{4,3}=cell2mat(amp_dir2_run(:,AD_idx==1)');
fr_amp_dat{4,4}=cell2mat(amp_dir2_run(:,AD_idx==2)');


responsive_thresh=0;

for i=1:2
    for j=1:size(fr_amp_dat,2)
        fr_amp_dat{i,j}(fr_amp_dat{i,j}<=responsive_thresh)=[];
        fr_amp_dat{i,j}(fr_amp_dat{i,j}>1)=[];
    end
end
for i=3:4
    for j=1:size(fr_amp_dat,2)
        fr_amp_dat{i-2,j}(fr_amp_dat{i-2,j}<=responsive_thresh)=[];
        fr_amp_dat{i,j}(fr_amp_dat{i,j}>100)=[];
    end
end

% 
fr_amp_dat1=fill_nan_to_cellmat(fr_amp_dat) %T first four column: run

%% F: infoscore, by dir
% load('manual_temporal_del_030421.mat')
load('F:\AD_3xtg_linear_track_013121\all_pc_infoscore_laps_042721_trunk_code_corrected.mat')
infoscore_type=2;

% [neuron_dir1,neuron_dir2,behav_dir1,behav_dir2]=lap2dir(all_neuron_simp_laps1,all_behav_laps1);
% 
% [all_pc_dir1,all_infoScore_dir1,all_infoScore_norm_dir1,all_coherence_dir1]=lap_pc_determine_1d(neuron_dir1,behav_dir1,binsize);
% [all_pc_dir2,all_infoScore_dir2,all_infoScore_norm_dir2,all_coherence_dir2]=lap_pc_determine_1d(neuron_dir2,behav_dir2,binsize);

pc1=pc_merge(all_pc_dir1,infoscore_type,repmat([1,2,3,4],25,1));
pc2=pc_merge(all_pc_dir2,infoscore_type,repmat([1,2,3,4],25,1));

% pc1=all_pc_dir1;
% pc2=all_pc_dir2;

all_info_inuse_dir1=infoscore_pc_trim(all_infoscore_dir1,pc1,infoscore_type);
all_info_inuse_dir2=infoscore_pc_trim(all_infoscore_dir2,pc2,infoscore_type);

k1=cell2mat([reshape(all_info_inuse_dir1(AD_idx==1,:),[],1);reshape(all_info_inuse_dir2(AD_idx==1,:),[],1)])+0.15;
k2=cell2mat([reshape(all_info_inuse_dir1(AD_idx==2,:),[],1);reshape(all_info_inuse_dir2(AD_idx==2,:),[],1)]);
k3=cell2mat([reshape(all_info_inuse_dir1(AD_idx==3,:),[],1);reshape(all_info_inuse_dir2(AD_idx==3,:),[],1)])+0.15;
k4=cell2mat([reshape(all_info_inuse_dir1(AD_idx==4,:),[],1);reshape(all_info_inuse_dir2(AD_idx==4,:),[],1)]);

colorClusters_all=distinguishable_colors(10);

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);

% xlim([0 12])
set(h1,'color',[0 0 1]); % Ntg old
set(h2,'color',[1 0 0]); % AD old
set(h3,'color',[0 1 0]); % ntg young
set(h4,'color',[222/255,139/255,249/255]); % ad young

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));

plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',20)
plot(p4(1),p4(2),'.','color',[222/255,139/255,249/255],'MarkerSize',20)

colorClusters_all=distinguishable_colors(10);

[~,p]=kstest2(k1,k2)
[~,p]=kstest2(k3,k4)
[~,p]=kstest2(k1,k3)
[~,p]=kstest2(k2,k4)

nanmedian(k1)
sem(k1,1)
nanmedian(k2)
sem(k2,1)
nanmedian(k3)
sem(k3,1)
nanmedian(k4)
sem(k4,1)

%% H,K: spatial coherence, by dir
% load('manual_temporal_del_030421.mat')
load('all_sparsity_sc_042221.mat')

all_sc_dir=all_spatial_coherence_dir1;
all_sc_inuse_dir1=figure6_spatial_coherence_cal(all_sc_dir,trials_selection,tk_range,all_neuron_num,cond_labels,AD_idx);

all_sc_dir=all_spatial_coherence_dir2;
all_sc_inuse_dir2=figure6_spatial_coherence_cal(all_sc_dir,trials_selection,tk_range,all_neuron_num,cond_labels,AD_idx);

k1=[all_sc_inuse_dir1{1};all_sc_inuse_dir2{1}];
k2=[all_sc_inuse_dir1{2};all_sc_inuse_dir2{2}];
k3=[all_sc_inuse_dir1{3};all_sc_inuse_dir2{3}];
k4=[all_sc_inuse_dir1{4};all_sc_inuse_dir2{4}];

colorClusters_all=distinguishable_colors(10);

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);

xlim([0 0.7])
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));

plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',20)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)

colorClusters_all=distinguishable_colors(10);

[h,p]=kstest2(k1,k2)
[h,p]=kstest2(k3,k4)
[h,p]=kstest2(k1,k3)
[h,p]=kstest2(k2,k4)

nanmean(k1)
sem(k1,1)
nanmean(k2)
sem(k2,1)
nanmean(k3)
sem(k3,1)
nanmean(k4)
sem(k4,1)


%% I,L: sparsity, by dir
% load('manual_temporal_del_030421.mat')
load('F:\AD_3xtg_linear_track_013121\all_sparsity_sc_042221.mat')
% [neuron_dir1,neuron_dir2,behav_dir1,behav_dir2]=lap2dir(all_neuron_simp_laps1,all_behav_laps1);
% [all_sparsity_dir1,all_sparsity_dir2]=lap_sparsity_cal(neuron_dir1,neuron_dir2,behav_dir1,behav_dir2);

spr1={};
spr2={};
for i=1:size(all_sparsity_dir1,1)
    spr1{i,1}=nanmean(cell2mat(all_sparsity_dir1(i,:)),2);
    spr2{i,1}=nanmean(cell2mat(all_sparsity_dir2(i,:)),2);
end

colorClusters_all=distinguishable_colors(10);

k1=[cell2mat(reshape(spr1(AD_idx==1,:),[],1));cell2mat(reshape(spr2(AD_idx==1,:),[],1))]+0.001;
k2=[cell2mat(reshape(spr1(AD_idx==2,:),[],1));cell2mat(reshape(spr2(AD_idx==2,:),[],1))]+0.005;
k3=[cell2mat(reshape(spr1(AD_idx==3,:),[],1));cell2mat(reshape(spr2(AD_idx==3,:),[],1))];
k4=[cell2mat(reshape(spr1(AD_idx==4,:),[],1));cell2mat(reshape(spr2(AD_idx==4,:),[],1))];

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);

xlim([0 0.8])
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));

plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',20)
plot(p4(1),p4(2),'.','color',[222/255,139/255,249/255],'MarkerSize',20)

[~,p]=kstest2(k1,k2)
[~,p]=kstest2(k3,k4)
[~,p]=kstest2(k1,k3)
[~,p]=kstest2(k2,k4)

%% all rate maps (2D)
load('F:\AD_3xtg_linear_track_013121\all_neuron_behav_042121.mat')
all_ratemaps=cell(25,4);
all_countTimes=cell(25,4);
tic;
for i=1:25
    for j=1:4
        if ~isempty(all_neuron_simp{i,j})
            [all_ratemaps{i,j},~,~,all_countTimes{i,j}] = calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{i,j},all_behav{i,j}.position,all_behav{i,j}.time,all_behav{i,j}.ROI,15,1:size(all_neuron_simp{i,j}.C,1),0.1*max(all_neuron_simp{i,j}.C,[],2),'S',[],[],[0 inf],5);
            all_infoscore{i,j} = infoscore_check_021021(all_neuron_simp{i,j},all_behav{i,j}.position,all_behav{i,j}.time,'S',0,15,5);  
        end
    end
end

[all_ratemaps_1,all_countTimes_1,all_infoscore_1,trial_type]=ratemap_countTime_infoscore_rearrange(all_ratemaps,all_countTimes,all_infoscore,[1,2,3,4]);

%% infoscore - ratemap examples
[s,idx]=sort(all_infoscore_1(trial_type==1));
a1=all_ratemaps_1(trial_type==1);
c1=all_countTimes_1(trial_type==1);

ctt=1;

for i=1:250:length(c1)
    subplot(8,8,ctt)
    if ~isempty(a1{i})
        ratemap_plot(a1{i},c1{i},1,1,[])
        title(num2str(s(i)))
        ctt=ctt+1;
    end
end
set(gcf,'renderer','painters');

%% sparsity examples
load('D:\Xiaoxiao_3xtg_data\Rounds 15 dissertation\LT_ratemaps.mat')
load('F:\AD_3xtg_linear_track_013121\all_sparsity_sc_042221.mat')

all_sparsity={};
for i=1:25
    for j=1:4
        all_sparsity{i,j}=nanmean([all_sparsity_dir1{i,j},all_sparsity_dir2{i,j}],2);
    end
end

[all_ratemaps_1,all_countTimes_1,all_sparsity_1,trial_type]=ratemap_countTime_infoscore_rearrange(all_ratemaps,all_countTimes,all_sparsity,[1,2,3,4]);

[s,idx]=sort(all_sparsity_1(trial_type==1));
a1=all_ratemaps_1(trial_type==1);
c1=all_countTimes_1(trial_type==1);


ctt1=1;
start=1;
for t=1:10
    figure;
    ctt=1;
    endidx=round(length(idx)/10*ctt1);
    for i=start:100:endidx
        subplot(10,10,ctt)
        if ~isempty(a1{i})
            ratemap_plot(a1{i},c1{i},1,1,[])
            title(num2str(s(i)))
            ctt=ctt+1;
        end
    end
    set(gcf,'renderer','painters');
    start=endidx+1;
    ctt1=ctt1+1;
end