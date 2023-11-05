%% new figure 3 panels: locomotion, use pearson, ensemble
load('F:\AD_square_circle_results_092320\exp_info_etgAD_031621.mat')
% load('D:\AD_square_circle_results_092320\manual_temporal_check_020521.mat')
del_ind=cell(1,30);
AD_idx=[4 4 4 4 3 3 3 3 4 4 3 3 4 4 2 2 2 1 1 1 2 2 1 2 1 1 1 2 2 1 1 1];
folderName=unique(destination);

bad_mice=[7 25]; % batch 5 972 is fake
tk_range=[1:length(folderName)];
tk_range(bad_mice)=[];

cond_labels={
    'Circle1';
    'Square1';
    'Square2';
    'Circle2';
    'Circle3';
    'Square3';
    'Square4';
    'Circle4';
    };
behavname={
    'circle1_behav.mat';
    'square1_behav.mat';
    'square2_behav.mat';
    'circle2_behav.mat';
    'circle3_behav.mat';
    'square3_behav.mat';
    'square4_behav.mat';
    'circle4_behav.mat';
    }
folderName=unique(destination);

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


% all_pc=cell(32,8);
% all_infoscore=cell(32,8);
% all_infoscore_norm=cell(32,8);
% all_coherence=cell(32,8);
% tic;
% parfor i=1:32*8
%     try
%         [place_cells,infoScore,infoScore_norm,coh] = permutingSpike_adapt_040821(all_neuron_simp_trial{i},all_behav{i}.position,all_behav{i}.time,'S',0,15,5,'all',0.5);  
% 
%         all_pc{i}=place_cells;
%         all_infoscore{i}=infoScore;
%         all_infoscore_norm{i}=infoScore_norm;
%         all_coherence{i}=coh;
%     catch
%     end
% end
% toc;

% load('D:\Xiaoxiao_3xtg_data\Rounds 14\all_neuron_behav_050321.mat')
% load('F:\AD_square_circle_results_092320\all_pc_infoscore_coherence_sparsity_051023.mat')

load('F:\AD_square_circle_results_092320\all_neuron_behav_032921.mat')
load('F:\AD_square_circle_results_092320\all_pc_infoscore_coherence_042921_goodFig5.mat')
infoscore_type=2;

%% all velo cal
all_velo={};
all_velo_raw={};
for tk=1:length(folderName)
    for j=1:length(behavname)
        
        all_velo{tk,j}=behav_velo_cal(all_behav{tk,j},[],'r');
        
        ntime=all_neuron_simp{tk,j}.time(1:2:end);
        ntime=resample(ntime,size(all_neuron_simp{tk,j}.C,2),length(ntime));

        all_velo{tk,j}=interp1(all_behav{tk,j}.time(1:end-1),all_velo{tk,j},ntime);
        all_velo_raw{tk,j}=all_velo{tk,j};
        h=fspecial('gaussian',[1,15],2); % 1 sec gauss
        all_velo{tk,j}=conv(all_velo{tk,j},h,'same');

    end
end

%% PANEL NEW: velo check for different types of mice

all_velo_test=[];
for tk=tk_range
    for j=trials_selection(tk,:)
        all_velo_test(tk,j)=nanmean(all_velo_raw{tk,j});
    end
end
all_velo_test(all_velo_test==0)=nan;

AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;

avm1=all_velo_test(AD_idx_t==1,:); % Ntg old
avm2=all_velo_test(AD_idx_t==2,:); % ad old
avm3=all_velo_test(AD_idx_t==3,:); % Ntg old
avm4=all_velo_test(AD_idx_t==4,:); % Ntg old

all_velo_test_type={reshape(avm3,size(avm3,1)*size(avm3,2),1),reshape(avm4,size(avm4,1)*size(avm4,2),1),reshape(avm1,size(avm1,1)*size(avm1,2),1),reshape(avm2,size(avm2,1)*size(avm2,2),1)};
all_velo_test_type_m=fill_nan_to_cellmat(all_velo_test_type)

all_velo_test_type_m{1}(all_velo_test_type_m{1}==0)=nan;
nanmax(all_velo_test_type_m{1}(:))

velo_thresh_by_type=nanmean(all_velo_test_type_m{1},1); % Ntg old, AD old, Ntg young, AD young

%% panel A: EXAMPLES
mice_idx=[4 11 23 19]+2; % Ntg young, AD young, Ntg old, AD old
thresh_selection=[3 4 2 1];
cond_idx=[ % circle, square
    4,2;
    1,2;
    1,2;
    1,2
    ];
ctt=1;
idx_col=2;
del_idx=del_ind;
neuron_nums=[];
for tk=mice_idx
    cd(folderName{tk});
      
    load('all_neurons.mat');
    nC=neurons{cond_idx(ctt,idx_col)}.C;
    nC(del_idx{tk},:)=[];
    neuron_nums(ctt)=size(nC,1);
    
    h1=subplot(3,4,ctt)
    imagesc(nC)
    colormap(h1,hot);
    caxis([0 70])
    
    h2=subplot(3,4,ctt+4)
    velo_all={};
    
    load(behavname{cond_idx(ctt,idx_col)});
    [velo]=behav_velo_cal(behav,[],'r');
    
    velo_all_m=velo(1:2:end);
    
%     h=ones(30,1)/30; % 1 sec avg
    h=fspecial('gaussian',[1,30],4); % 2 sec gauss
    velo_all_m=conv(velo_all_m,h,'same');
    velo_all_m=velo_all_m>=velo_thresh_by_type(thresh_selection(ctt)); %10mm/sec
    
    imagesc(~velo_all_m');
    colormap(h2,'gray');
    
    h3=subplot(3,4,ctt+8)
    h=fspecial('gaussian',[1,150],round((150-1)/4)); % 2 sec gauss
    velo_all_m1=conv(velo(1:2:end),h,'same');
    plot(zscore(velo_all_m1));
    hold on;
    sumC=sum(nC,1);
    sumC=conv(sumC,h,'same');
    plot(zscore(sumC));
    
    ctt=ctt+1;
end
set(gcf,'renderer','painters');


%% 3: amp-activate rate (ratio of running frames in nearby 150 frames) corr, pearson this time
% load('D:\AD_square_circle_results_092320\manual_temporal_check_add_lowfr.mat')
h_gauss=fspecial('gaussian',[1,150],round((150-1)/4));
h_mean=fspecial('average',[1,150]);

corr_fr_runEpoch_all={};
corr_fr_runEpoch_all_m={};
tic;
for tk=1:length(foldername)
    for j=trials_selection(tk,:)
%         nC=all_neuron_simp{tk,j}.C;
        nC=all_neuron_simp{tk,j}.C;

        av_1=fillmissing(all_velo{tk,j},'linear');
        av_1=conv(av_1,h_gauss,'same');
        
        corr_fr_runEpoch=[];
        
        for k=1:size(nC,1)
            fr_smooth1=conv(nC(k,:),h_gauss,'same');
            [corrt]=corrcoef(zscore(fr_smooth1'),zscore(av_1));
            corr_fr_runEpoch(k)=corrt(2);
        end
        
        corr_fr_runEpoch_all{tk,j}=corr_fr_runEpoch';
%         lag_fr_runEpoch_all{tk,j}=lag_fr_runEpoch';
%         corr_fr_nrunEpoch_all{tk,j}=corr_fr_nrunEpoch';
        fr_smooth3=conv(sum(nC,1),h_gauss,'same');
%         fr_smooth2=C_to_peakS_with_thresh(fr_smooth3,0.1*nanmax(fr_smooth3));
%         fr_smooth2=movmean(fr_smooth2>0,151);
%         fr_smooth2=conv(fr_smooth2,h_gauss,'same');
        fr_smooth2=fr_smooth3;
        [corrt]=corrcoef(zscore(fr_smooth2'),zscore(av_1),'rows','complete');
        corr_fr_runEpoch_all_m{tk,j}=corrt(2);
    end
    toc;
end

corr_fr_both_all={};
corr_fr_both_m={};

for tk=1:length(foldername)
    trial_select_cir=[];
    trial_select_sqr=[];
    for j=1:length(trials_selection(tk,:))
        if ~isempty(strfind(cond_labels{trials_selection(tk,j)},'Circle'))
            trial_select_cir=[trial_select_cir,trials_selection(tk,j)];
        else
            trial_select_sqr=[trial_select_sqr,trials_selection(tk,j)];
        end
    end
    corr_fr_both_all{tk,1}=[nanmean(cell2mat(corr_fr_runEpoch_all(tk,trial_select_cir)),2);nanmean(cell2mat(corr_fr_runEpoch_all(tk,trial_select_sqr)),2)];

    corr_fr_both_m{tk,1}=nanmean([nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_cir)));nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_sqr)))]);
end

corr_fr_both_m2_amp={cell2mat(corr_fr_both_all(AD_idx==3,:)),cell2mat(corr_fr_both_all(AD_idx==4,:)),cell2mat(corr_fr_both_all(AD_idx==1,:)),cell2mat(corr_fr_both_all(AD_idx==2,:))};
corr_fr_both_m2_amp=fill_nan_to_cellmat(corr_fr_both_m2_amp);

corr_fr_both_mt_amp={cell2mat(corr_fr_both_m(AD_idx==3)),cell2mat(corr_fr_both_m(AD_idx==4)),cell2mat(corr_fr_both_m(AD_idx==1)),cell2mat(corr_fr_both_m(AD_idx==2))};
corr_fr_both_mt1_amp=fill_nan_to_cellmat(corr_fr_both_mt_amp);

k1=corr_fr_both_m2_amp{1}(:,3);
k2=corr_fr_both_m2_amp{1}(:,4);
k3=corr_fr_both_m2_amp{1}(:,1);
k4=corr_fr_both_m2_amp{1}(:,2);

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
% xlim([0 0.8])
set(h1,'color',colorClusters_all(1,:));
set(h2,'color',colorClusters_all(2,:));
set(h3,'color',colorClusters_all(3,:));
set(h4,'color',[222/255,139/255,249/255]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));

plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',20)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)

[h,p]=kstest2(k1,k2)
[h,p]=kstest2(k3,k4)
[h,p]=kstest2(k1,k3)
[h,p]=kstest2(k2,k4)

ranksum(k1,k2)
ranksum(k3,k4)
ranksum(k1,k3)
ranksum(k2,k4)

% special: per mice
AD_idx_t=AD_idx;AD_idx_t(bad_mice)=-1;
corr_fr_both_m=cell2mat(corr_fr_both_m);
dat=fill_nan_to_cellmat({corr_fr_both_m(AD_idx_t==1),corr_fr_both_m(AD_idx_t==2),corr_fr_both_m(AD_idx_t==3),corr_fr_both_m(AD_idx_t==4)});
%% 5: amp-activate rate (ratio of running frames in nearby 150 frames) corr, place cells, pearson this time
% load('D:\AD_square_circle_results_092320\manual_temporal_check_add_lowfr.mat')
h_gauss=fspecial('gaussian',[1,150],round((150-1)/4));
h_mean=fspecial('average',[1,150]);

corr_fr_runEpoch_all={};
corr_fr_runEpoch_all_m={};
lag_fr_runEpoch_all={};
lag_fr_runEpoch_all_m={};
tic;
infoscore_type=2;
for tk=tk_range
    for j=trials_selection(tk,:)
%         nC=all_neuron_simp{tk,j}.C;
        nC=all_neuron_simp{tk,j}.C(all_pc{tk,j}{infoscore_type},:);

        av_1=fillmissing(all_velo{tk,j},'linear');
        av_1=conv(av_1,h_gauss,'same');
        
        corr_fr_runEpoch=[];      
        for k=1:size(nC,1)
            fr_smooth1=conv(nC(k,:),h_gauss,'same');
            [corrt]=corrcoef(zscore(fr_smooth1'),zscore(av_1));
            corr_fr_runEpoch(k)=corrt(2);

        end
        
        corr_fr_runEpoch_all{tk,j}=corr_fr_runEpoch';
        fr_smooth3=conv(sum(nC,1),h_gauss,'same');
        fr_smooth2=fr_smooth3;
        [corrt]=corrcoef(zscore(fr_smooth2'),zscore(av_1));
        corr_fr_runEpoch_all_m{tk,j}=corrt(2);
    end
    toc;
end


corr_fr_both_all={};
corr_fr_both_m={};

for tk=tk_range
    trial_select_cir=[];
    trial_select_sqr=[];
    for j=1:length(trials_selection(tk,:))
        if ~isempty(strfind(cond_labels{trials_selection(tk,j)},'Circle'))
            trial_select_cir=[trial_select_cir,trials_selection(tk,j)];
        else
            trial_select_sqr=[trial_select_sqr,trials_selection(tk,j)];
        end
    end
    corr_fr_both_all{tk,1}=[cell2mat(corr_fr_runEpoch_all(tk,trials_selection(tk,:))')];
    corr_fr_both_m{tk,1}=nanmean([nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_cir)));nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_sqr)))]);
end

corr_fr_both_m2_amp_pc={cell2mat(corr_fr_both_all(AD_idx==3,:)),cell2mat(corr_fr_both_all(AD_idx==4,:)),cell2mat(corr_fr_both_all(AD_idx==1,:)),cell2mat(corr_fr_both_all(AD_idx==2,:))};
corr_fr_both_m2_amp_pc=fill_nan_to_cellmat(corr_fr_both_m2_amp_pc);

corr_fr_both_mt_amp_pc={cell2mat(corr_fr_both_m(AD_idx==3)),cell2mat(corr_fr_both_m(AD_idx==4)),cell2mat(corr_fr_both_m(AD_idx==1)),cell2mat(corr_fr_both_m(AD_idx==2))};
corr_fr_both_mt1_amp_pc=fill_nan_to_cellmat(corr_fr_both_mt_amp_pc);

k1=corr_fr_both_m2_amp_pc{1}(:,3);
k2=corr_fr_both_m2_amp_pc{1}(:,4);
k3=corr_fr_both_m2_amp_pc{1}(:,1);
k4=corr_fr_both_m2_amp_pc{1}(:,2);

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
% xlim([0 0.8])
set(h1,'color',colorClusters_all(1,:));
set(h2,'color',colorClusters_all(2,:));
set(h3,'color',colorClusters_all(3,:));
set(h4,'color',[222/255,139/255,249/255]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));

plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',20)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)

[h,p]=kstest2(k1,k2)
[h,p]=kstest2(k3,k4)
[h,p]=kstest2(k1,k3)
[h,p]=kstest2(k2,k4)

ranksum(k1,k2)
ranksum(k3,k4)
ranksum(k1,k3)
ranksum(k2,k4)

nanmedian(k3)
nanmedian(k4)
nanmedian(k1)
nanmedian(k2)
