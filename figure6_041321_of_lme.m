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

all_pc=cell(32,8);
all_infoscore=cell(32,8);
all_infoscore_norm=cell(32,8);
all_coherence=cell(32,8);
tic;
parfor i=1:32*8
    try
        [place_cells,infoScore,infoScore_norm,coh] = permutingSpike_adapt_040821(all_neuron_simp_trial{i},all_behav{i}.position,all_behav{i}.time,'S',0,10,5,'all',0.3);  

        all_pc{i}=place_cells;
        all_infoscore{i}=infoScore;
        all_infoscore_norm{i}=infoScore_norm;
        all_coherence{i}=coh;
    catch
    end
end
toc;

infoscore_type=2;

%% all velo cal
all_velo={};
all_velo_raw={};
for tk=1:length(folderName)
    load([folderName{tk},'\','neuronIndividuals_new_trial.mat']);
    for j=1:length(behavname)
        load([folderName{tk},'\',behavname{j}]);
        
        all_velo{tk,j}=behav_velo_cal(behav,[],'r');
        
        ntime=neuronIndividuals_new{j}.time(1:2:end);
        ntime=resample(ntime,size(neuronIndividuals_new{j}.C,2),length(ntime));

        all_velo_raw{tk,j}=all_velo{tk,j};
        all_velo{tk,j}=interp1(behav.time(1:end-1),all_velo{tk,j},ntime);
        
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
for tk=mice_idx
    cd(folderName{tk});
      
    load('all_neurons.mat');
    nC=neurons{cond_idx(ctt,idx_col)}.C;
    nC(del_idx{tk},:)=[];
    
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

%% 3: fr-activate rate (ratio of running frames in nearby 150 frames) corr, pearson this time
load('F:\AD_square_circle_results_092320\velocity_033021.mat');
load('F:\AD_square_circle_results_092320\all_neuron_behav_trial_040321.mat')
% load('D:\AD_square_circle_results_092320\manual_temporal_check_add_lowfr.mat')
h_gauss=fspecial('gaussian',[1,150],round((150-1)/4));
h_mean=fspecial('average',[1,150]);

corr_fr_runEpoch_all={};
corr_fr_runEpoch_all_m={};
lag_fr_runEpoch_all={};
lag_fr_runEpoch_all_m={};
tic;
for tk=tk_range
    for j=trials_selection(tk,:)
%         nC=all_neuron_simp{tk,j}.C;
        nS=all_neuron_simp_trial{tk,j}.S(all_pc{tk,j}{infoscore_type},:);

        av_1=fillmissing(all_velo{tk,j},'linear');
        av_1=conv(av_1,h_gauss,'same');
        
        corr_fr_runEpoch=[];
        lag_fr_runEpoch=[];
        
        for k=1:size(nS,1)
            frS=conv(nS(k,:)>0,h_mean,'same');
            fr_smooth1=conv(frS,h_gauss,'same');
            [corrt]=corrcoef(zscore(fr_smooth1'),zscore(av_1));
            corr_fr_runEpoch(k)=corrt(2);
            if corr_fr_runEpoch(k)>0
                lag_fr_runEpoch(k)=0;
            else
                lag_fr_runEpoch(k)=0;
            end
        end
        
        corr_fr_runEpoch_all{tk,j}=corr_fr_runEpoch';
%         lag_fr_runEpoch_all{tk,j}=lag_fr_runEpoch';
%         corr_fr_nrunEpoch_all{tk,j}=corr_fr_nrunEpoch';
        frS3=conv(double(sum(nS,1)>0),h_mean,'same');
        fr_smooth3=conv(frS3,h_gauss,'same');
%         fr_smooth2=C_to_peakS_with_thresh(fr_smooth3,0.1*nanmax(fr_smooth3));
%         fr_smooth2=movmean(fr_smooth2>0,151);
%         fr_smooth2=conv(fr_smooth2,h_gauss,'same');
        fr_smooth2=fr_smooth3;
        [corrt]=corrcoef(zscore(fr_smooth2'),zscore(av_1));
        corr_fr_runEpoch_all_m{tk,j}=corrt(2);
    end
    toc;
end

corr_fr_cir={};
corr_fr_sqr={};
corr_fr_both={};
corr_fr_both_all={};
corr_fr_both_cellmean={};
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
%     corr_fr_both_m{tk,1}=[nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_cir)));nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_sqr)))]
%     corr_fr_both_m{tk,1}=[cell2mat(corr_fr_runEpoch_all_m(tk,trials_selection(tk,:))')];
end

corr_fr_both_m2_fr={cell2mat(corr_fr_both_all(AD_idx==3,:)),cell2mat(corr_fr_both_all(AD_idx==4,:)),cell2mat(corr_fr_both_all(AD_idx==1,:)),cell2mat(corr_fr_both_all(AD_idx==2,:))};
corr_fr_both_m2_fr=fill_nan_to_cellmat(corr_fr_both_m2_fr);

corr_fr_both_mt_fr={cell2mat(corr_fr_both_m(AD_idx==3)),cell2mat(corr_fr_both_m(AD_idx==4)),cell2mat(corr_fr_both_m(AD_idx==1)),cell2mat(corr_fr_both_m(AD_idx==2))};
corr_fr_both_mt1_fr=fill_nan_to_cellmat(corr_fr_both_mt_fr);

%% 4: amp-activate rate (ratio of running frames in nearby 150 frames) corr, pearson this time
load('F:\AD_square_circle_results_092320\velocity_033021.mat');
% load('F:\AD_square_circle_results_092320\all_neuron_behav_trial_040321.mat')
% load('D:\AD_square_circle_results_092320\manual_temporal_check_add_lowfr.mat')
h_gauss=fspecial('gaussian',[1,150],round((150-1)/4));
h_mean=fspecial('average',[1,150]);

corr_fr_runEpoch_all={};
corr_fr_runEpoch_all_m={};
lag_fr_runEpoch_all={};
lag_fr_runEpoch_all_m={};
tic;
for tk=tk_range
    for j=trials_selection(tk,:)
%         nC=all_neuron_simp{tk,j}.C;
        nC=all_neuron_simp_trial{tk,j}.C(all_pc{tk,j}{infoscore_type},:);

        av_1=fillmissing(all_velo{tk,j},'linear');
        av_1=conv(av_1,h_gauss,'same');
        
        corr_fr_runEpoch=[];
        lag_fr_runEpoch=[];
        
        for k=1:size(nC,1)
            fr_smooth1=conv(nC(k,:),h_gauss,'same');
            [corrt]=corrcoef(zscore(fr_smooth1'),zscore(av_1));
            corr_fr_runEpoch(k)=corrt(2);
            if corr_fr_runEpoch(k)>0
                lag_fr_runEpoch(k)=0;
            else
                lag_fr_runEpoch(k)=0;
            end
        end
        
        corr_fr_runEpoch_all{tk,j}=corr_fr_runEpoch';
%         lag_fr_runEpoch_all{tk,j}=lag_fr_runEpoch';
%         corr_fr_nrunEpoch_all{tk,j}=corr_fr_nrunEpoch';
        fr_smooth3=conv(sum(nC,1),h_gauss,'same');
%         fr_smooth2=C_to_peakS_with_thresh(fr_smooth3,0.1*nanmax(fr_smooth3));
%         fr_smooth2=movmean(fr_smooth2>0,151);
%         fr_smooth2=conv(fr_smooth2,h_gauss,'same');
        fr_smooth2=fr_smooth3;
        [corrt]=corrcoef(zscore(fr_smooth2'),zscore(av_1));
        corr_fr_runEpoch_all_m{tk,j}=corrt(2);
    end
    toc;
end

corr_fr_cir={};
corr_fr_sqr={};
corr_fr_both={};
corr_fr_both_all={};
corr_fr_both_cellmean={};
corr_fr_both_m={};

type_idx={};
trial_idx={};
m_idx={};
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
%     corr_fr_both_m{tk,1}=[nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_cir)));nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_sqr)))]
%     corr_fr_both_m{tk,1}=[cell2mat(corr_fr_runEpoch_all_m(tk,trials_selection(tk,:))')];
    type_idx{tk,1}=ones(length(corr_fr_both_all{tk,1}),1)*AD_idx(tk);
    trial_idx{tk,1}=[ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,1)}),1)*1;ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,2)}),1)*2;ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,3)}),1)*3;ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,4)}),1)*4]; 
    m_idx{tk,1}=ones(length(corr_fr_both_all{tk,1}),1)*tk;


end

corr_fr_both_m2_amp={cell2mat(corr_fr_both_all(AD_idx==3,:)),cell2mat(corr_fr_both_all(AD_idx==4,:)),cell2mat(corr_fr_both_all(AD_idx==1,:)),cell2mat(corr_fr_both_all(AD_idx==2,:))};
corr_fr_both_m2_amp=fill_nan_to_cellmat(corr_fr_both_m2_amp);

corr_fr_both_mt_amp={cell2mat(corr_fr_both_m(AD_idx==3)),cell2mat(corr_fr_both_m(AD_idx==4)),cell2mat(corr_fr_both_m(AD_idx==1)),cell2mat(corr_fr_both_m(AD_idx==2))};
corr_fr_both_mt1_amp=fill_nan_to_cellmat(corr_fr_both_mt_amp);

var=cell2mat(corr_fr_both_all);
type_idx1=cell2mat(type_idx);
trial_idx1=cell2mat(trial_idx);
m_idx1=cell2mat(m_idx);

idx_nan=isnan(var);
var(idx_nan)=[];
type_idx1(idx_nan)=[];
trial_idx1(idx_nan)=[];
m_idx1(idx_nan)=[];

type_idx11=type_idx1;
type_idx11(type_idx1==1)=0; % reference group is 1 or 0
[lme1,anova_lme]=calculate_lme_011221(var,{categorical(type_idx11)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

type_idx12=type_idx1;
type_idx12(type_idx1==2)=0; % reference group is 1 or 0
[lme2,anova_lme]=calculate_lme_011221(var,{categorical(type_idx12)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

type_idx13=type_idx1;
type_idx13(type_idx1==3)=0; % reference group is 1 or 0
[lme3,anova_lme]=calculate_lme_011221(var,{categorical(type_idx13)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

type_idx14=type_idx1;
type_idx14(type_idx1==4)=0; % reference group is 1 or 0
[lme4,anova_lme]=calculate_lme_011221(var,{categorical(type_idx14)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

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
        lag_fr_runEpoch=[];
        
        for k=1:size(nC,1)
            fr_smooth1=conv(nC(k,:),h_gauss,'same');
            [corrt]=corrcoef(zscore(fr_smooth1'),zscore(av_1));
            corr_fr_runEpoch(k)=corrt(2);
            if corr_fr_runEpoch(k)>0
                lag_fr_runEpoch(k)=0;
            else
                lag_fr_runEpoch(k)=0;
            end
        end
        
        corr_fr_runEpoch_all{tk,j}=corr_fr_runEpoch';
%         lag_fr_runEpoch_all{tk,j}=lag_fr_runEpoch';
%         corr_fr_nrunEpoch_all{tk,j}=corr_fr_nrunEpoch';
        fr_smooth3=conv(sum(nC,1),h_gauss,'same');
%         fr_smooth2=C_to_peakS_with_thresh(fr_smooth3,0.1*nanmax(fr_smooth3));
%         fr_smooth2=movmean(fr_smooth2>0,151);
%         fr_smooth2=conv(fr_smooth2,h_gauss,'same');
        fr_smooth2=fr_smooth3;
        [corrt]=corrcoef(zscore(fr_smooth2'),zscore(av_1));
        corr_fr_runEpoch_all_m{tk,j}=corrt(2);
    end
    toc;
end

corr_fr_cir={};
corr_fr_sqr={};
corr_fr_both={};
corr_fr_both_all={};
corr_fr_both_cellmean={};
corr_fr_both_m={};

type_idx={};
trial_idx={};
m_idx={};
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
%     corr_fr_both_m{tk,1}=[nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_cir)));nanmean(cell2mat(corr_fr_runEpoch_all_m(tk,trial_select_sqr)))]
%     corr_fr_both_m{tk,1}=[cell2mat(corr_fr_runEpoch_all_m(tk,trials_selection(tk,:))')];
    type_idx{tk,1}=ones(length(corr_fr_both_all{tk,1}),1)*AD_idx(tk);
    trial_idx{tk,1}=[ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,1)}),1)*1;ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,2)}),1)*2;ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,3)}),1)*3;ones(length(corr_fr_runEpoch_all{tk,trials_selection(tk,4)}),1)*4]; 
    m_idx{tk,1}=ones(length(corr_fr_both_all{tk,1}),1)*tk;


end

corr_fr_both_m2_amp={cell2mat(corr_fr_both_all(AD_idx==3,:)),cell2mat(corr_fr_both_all(AD_idx==4,:)),cell2mat(corr_fr_both_all(AD_idx==1,:)),cell2mat(corr_fr_both_all(AD_idx==2,:))};
corr_fr_both_m2_amp=fill_nan_to_cellmat(corr_fr_both_m2_amp);

corr_fr_both_mt_amp={cell2mat(corr_fr_both_m(AD_idx==3)),cell2mat(corr_fr_both_m(AD_idx==4)),cell2mat(corr_fr_both_m(AD_idx==1)),cell2mat(corr_fr_both_m(AD_idx==2))};
corr_fr_both_mt1_amp=fill_nan_to_cellmat(corr_fr_both_mt_amp);

var=cell2mat(corr_fr_both_all);
type_idx1=cell2mat(type_idx);
trial_idx1=cell2mat(trial_idx);
m_idx1=cell2mat(m_idx);

idx_nan=isnan(var);
var(idx_nan)=[];
type_idx1(idx_nan)=[];
trial_idx1(idx_nan)=[];
m_idx1(idx_nan)=[];

type_idx11=type_idx1;
type_idx11(type_idx1==1)=0; % reference group is 1 or 0
[lme1,anova_lme]=calculate_lme_011221(var,{categorical(type_idx11)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

type_idx12=type_idx1;
type_idx12(type_idx1==2)=0; % reference group is 1 or 0
[lme2,anova_lme]=calculate_lme_011221(var,{categorical(type_idx12)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

type_idx13=type_idx1;
type_idx13(type_idx1==3)=0; % reference group is 1 or 0
[lme3,anova_lme]=calculate_lme_011221(var,{categorical(type_idx13)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');

type_idx14=type_idx1;
type_idx14(type_idx1==4)=0; % reference group is 1 or 0
[lme4,anova_lme]=calculate_lme_011221(var,{categorical(type_idx14)},{categorical(m_idx1)},'y~fix1+(fix1|random1)');
