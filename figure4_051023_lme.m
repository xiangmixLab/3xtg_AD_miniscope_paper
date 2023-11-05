%% figure 4 linear track
load('all_neuron_behav_042121.mat')
load('all_neuron_behav_laps_042221.mat')
load('all_sparsity_sc_042221.mat')
load('all_pc_infoscore_laps_042221.mat')
%% new E, F: fr, amp, both dir merge
% load('all_neuron_behav_waterrm_041621.mat');
load('F:\AD_3xtg_linear_track_013121\dat_info.mat');

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

fr_dir1_run={};
fr_dir2_run={};

amp_dir1_run={};
amp_dir2_run={};

del_idx=cell(25,1);
[neuron_dir]=laps_avg_neuroDat(all_neuron_simp_laps1,all_behav_laps1);

for tk=tk_range
    
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

%% E,F: fr lme
k1=[fr_amp_dat1{1,1}(:,1);fr_amp_dat1{2,1}(:,1)];
k2=[fr_amp_dat1{1,1}(:,2);fr_amp_dat1{2,1}(:,2)];
k3=[fr_amp_dat1{1,1}(:,3);fr_amp_dat1{2,1}(:,3)];
k4=[fr_amp_dat1{1,1}(:,4);fr_amp_dat1{2,1}(:,4)];

ranksum(k1,k2)
ranksum(k3,k4)
ranksum(k1,k3)
ranksum(k2,k4)

type_dir1_idx={};
trial_dir1_idx={};
m_dir1_idx={};

type_dir2_idx={};
trial_dir2_idx={};
m_dir2_idx={};

for i=tk_range
    type_dir1_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*AD_idx(i);
    trial_dir1_idx{i,1}=ones(size(fr_dir1_run{i},1),1);
    m_dir1_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*i;
    
    type_dir2_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*AD_idx(i);
    trial_dir2_idx{i,1}=ones(size(fr_dir1_run{i},1),1);
    m_dir2_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*i;
end

trial_dir1_idx_all={cell2mat(type_dir1_idx(AD_idx==1)),cell2mat(type_dir1_idx(AD_idx==2)),cell2mat(type_dir1_idx(AD_idx==3)),cell2mat(type_dir1_idx(AD_idx==4))};
trial_dir2_idx_all={cell2mat(type_dir2_idx(AD_idx==1)),cell2mat(type_dir2_idx(AD_idx==2)),cell2mat(type_dir2_idx(AD_idx==3)),cell2mat(type_dir2_idx(AD_idx==4))};

type_dir1_idx_all={cell2mat(trial_dir1_idx(AD_idx==1)),cell2mat(trial_dir1_idx(AD_idx==2)),cell2mat(trial_dir1_idx(AD_idx==3)),cell2mat(trial_dir1_idx(AD_idx==4))};
type_dir2_idx_all={cell2mat(trial_dir2_idx(AD_idx==1)),cell2mat(trial_dir2_idx(AD_idx==2)),cell2mat(trial_dir2_idx(AD_idx==3)),cell2mat(trial_dir2_idx(AD_idx==4))};

m_dir1_idx_all={cell2mat(m_dir1_idx(AD_idx==1)),cell2mat(m_dir1_idx(AD_idx==2)),cell2mat(m_dir1_idx(AD_idx==3)),cell2mat(m_dir1_idx(AD_idx==4))};
m_dir2_idx_all={cell2mat(m_dir2_idx(AD_idx==1)),cell2mat(m_dir2_idx(AD_idx==2)),cell2mat(m_dir2_idx(AD_idx==3)),cell2mat(m_dir2_idx(AD_idx==4))};

trial_dir1_idx_all_m=fill_nan_to_cellmat(trial_dir1_idx_all);
trial_dir2_idx_all_m=fill_nan_to_cellmat(trial_dir2_idx_all);
type_dir1_idx_all_m=fill_nan_to_cellmat(type_dir1_idx_all);
type_dir2_idx_all_m=fill_nan_to_cellmat(type_dir2_idx_all);
m_dir1_idx_all_m=fill_nan_to_cellmat(m_dir1_idx_all);
m_dir2_idx_all_m=fill_nan_to_cellmat(m_dir2_idx_all);

var=[k1;k2;k3;k4];
% trial_idx_m=[[trial_dir1_idx_all_m{1}(:,1);trial_dir2_idx_all_m{1}(:,1)];[trial_dir1_idx_all_m{1}(:,2);trial_dir2_idx_all_m{1}(:,2)];[trial_dir1_idx_all_m{1}(:,3);trial_dir2_idx_all_m{1}(:,3)];[trial_dir1_idx_all_m{1}(:,4);trial_dir2_idx_all_m{1}(:,4)]];
% type_idx_m=[[type_dir1_idx_all_m{1}(:,1);type_dir2_idx_all_m{1}(:,1)];[type_dir1_idx_all_m{1}(:,2);type_dir2_idx_all_m{1}(:,2)];[type_dir1_idx_all_m{1}(:,3);type_dir2_idx_all_m{1}(:,3)];[type_dir1_idx_all_m{1}(:,4);type_dir2_idx_all_m{1}(:,4)]];
% m_idx_m=[[m_dir1_idx_all_m{1}(:,1);m_dir2_idx_all_m{1}(:,1)];[m_dir1_idx_all_m{1}(:,2);m_dir2_idx_all_m{1}(:,2)];[m_dir1_idx_all_m{1}(:,3);m_dir2_idx_all_m{1}(:,3)];[m_dir1_idx_all_m{1}(:,4);m_dir2_idx_all_m{1}(:,4)]];
trial_idx_m=[ones(length(k1)/2,1);ones(length(k1)/2,1)*2;ones(length(k2)/2,1);ones(length(k2)/2,1)*2;ones(length(k3)/2,1);ones(length(k3)/2,1)*2;ones(length(k4)/2,1);ones(length(k4)/2,1)*2];
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[[m_dir1_idx_all_m{1}(:,1);m_dir2_idx_all_m{1}(:,1)];[m_dir1_idx_all_m{1}(:,2);m_dir2_idx_all_m{1}(:,2)];[m_dir1_idx_all_m{1}(:,3);m_dir2_idx_all_m{1}(:,3)];[m_dir1_idx_all_m{1}(:,4);m_dir2_idx_all_m{1}(:,4)]];

nanvar_idx=isnan(var);
var(nanvar_idx)=[];
trial_idx_m(nanvar_idx)=[];
type_idx_m(nanvar_idx)=[];
m_idx_m(nanvar_idx)=[];

type_idx_m1=type_idx_m; % ntg old
type_idx_m1(type_idx_m1==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

%% E,F: AMP LME
k1=[fr_amp_dat1{3,1}(:,1);fr_amp_dat1{4,1}(:,1)];
k2=[fr_amp_dat1{3,1}(:,2);fr_amp_dat1{4,1}(:,2)];
k3=[fr_amp_dat1{3,1}(:,3);fr_amp_dat1{4,1}(:,3)];
k4=[fr_amp_dat1{3,1}(:,4);fr_amp_dat1{4,1}(:,4)];

ranksum(k1,k2)
ranksum(k3,k4)
ranksum(k1,k3)
ranksum(k2,k4)

type_dir1_idx={};
trial_dir1_idx={};
m_dir1_idx={};

type_dir2_idx={};
trial_dir2_idx={};
m_dir2_idx={};

for i=tk_range
    type_dir1_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*AD_idx(i);
    trial_dir1_idx{i,1}=ones(size(fr_dir1_run{i},1),1);
    m_dir1_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*i;
    
    type_dir2_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*AD_idx(i);
    trial_dir2_idx{i,1}=ones(size(fr_dir1_run{i},1),1);
    m_dir2_idx{i,1}=ones(size(fr_dir1_run{i},1),1)*i;
end

trial_dir1_idx_all={cell2mat(type_dir1_idx(AD_idx==1)),cell2mat(type_dir1_idx(AD_idx==2)),cell2mat(type_dir1_idx(AD_idx==3)),cell2mat(type_dir1_idx(AD_idx==4))};
trial_dir2_idx_all={cell2mat(type_dir2_idx(AD_idx==1)),cell2mat(type_dir2_idx(AD_idx==2)),cell2mat(type_dir2_idx(AD_idx==3)),cell2mat(type_dir2_idx(AD_idx==4))};

type_dir1_idx_all={cell2mat(trial_dir1_idx(AD_idx==1)),cell2mat(trial_dir1_idx(AD_idx==2)),cell2mat(trial_dir1_idx(AD_idx==3)),cell2mat(trial_dir1_idx(AD_idx==4))};
type_dir2_idx_all={cell2mat(trial_dir2_idx(AD_idx==1)),cell2mat(trial_dir2_idx(AD_idx==2)),cell2mat(trial_dir2_idx(AD_idx==3)),cell2mat(trial_dir2_idx(AD_idx==4))};

m_dir1_idx_all={cell2mat(m_dir1_idx(AD_idx==1)),cell2mat(m_dir1_idx(AD_idx==2)),cell2mat(m_dir1_idx(AD_idx==3)),cell2mat(m_dir1_idx(AD_idx==4))};
m_dir2_idx_all={cell2mat(m_dir2_idx(AD_idx==1)),cell2mat(m_dir2_idx(AD_idx==2)),cell2mat(m_dir2_idx(AD_idx==3)),cell2mat(m_dir2_idx(AD_idx==4))};

trial_dir1_idx_all_m=fill_nan_to_cellmat(trial_dir1_idx_all);
trial_dir2_idx_all_m=fill_nan_to_cellmat(trial_dir2_idx_all);
type_dir1_idx_all_m=fill_nan_to_cellmat(type_dir1_idx_all);
type_dir2_idx_all_m=fill_nan_to_cellmat(type_dir2_idx_all);
m_dir1_idx_all_m=fill_nan_to_cellmat(m_dir1_idx_all);
m_dir2_idx_all_m=fill_nan_to_cellmat(m_dir2_idx_all);

var=[k1;k2;k3;k4];
% trial_idx_m=[[trial_dir1_idx_all_m{1}(:,1);trial_dir2_idx_all_m{1}(:,1)];[trial_dir1_idx_all_m{1}(:,2);trial_dir2_idx_all_m{1}(:,2)];[trial_dir1_idx_all_m{1}(:,3);trial_dir2_idx_all_m{1}(:,3)];[trial_dir1_idx_all_m{1}(:,4);trial_dir2_idx_all_m{1}(:,4)]];
% type_idx_m=[[type_dir1_idx_all_m{1}(:,1);type_dir2_idx_all_m{1}(:,1)];[type_dir1_idx_all_m{1}(:,2);type_dir2_idx_all_m{1}(:,2)];[type_dir1_idx_all_m{1}(:,3);type_dir2_idx_all_m{1}(:,3)];[type_dir1_idx_all_m{1}(:,4);type_dir2_idx_all_m{1}(:,4)]];
% m_idx_m=[[m_dir1_idx_all_m{1}(:,1);m_dir2_idx_all_m{1}(:,1)];[m_dir1_idx_all_m{1}(:,2);m_dir2_idx_all_m{1}(:,2)];[m_dir1_idx_all_m{1}(:,3);m_dir2_idx_all_m{1}(:,3)];[m_dir1_idx_all_m{1}(:,4);m_dir2_idx_all_m{1}(:,4)]];
trial_idx_m=[ones(length(k1)/2,1);ones(length(k1)/2,1)*2;ones(length(k2)/2,1);ones(length(k2)/2,1)*2;ones(length(k3)/2,1);ones(length(k3)/2,1)*2;ones(length(k4)/2,1);ones(length(k4)/2,1)*2];
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[[m_dir1_idx_all_m{1}(:,1);m_dir2_idx_all_m{1}(:,1)];[m_dir1_idx_all_m{1}(:,2);m_dir2_idx_all_m{1}(:,2)];[m_dir1_idx_all_m{1}(:,3);m_dir2_idx_all_m{1}(:,3)];[m_dir1_idx_all_m{1}(:,4);m_dir2_idx_all_m{1}(:,4)]];

nanvar_idx=isnan(var);
var(nanvar_idx)=[];
trial_idx_m(nanvar_idx)=[];
type_idx_m(nanvar_idx)=[];
m_idx_m(nanvar_idx)=[];

type_idx_m1=type_idx_m; % ntg old
type_idx_m1(type_idx_m1==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

%% infoscore LME
% suppose the k1-k4 as already been calculated

[var,trial_idx_m,type_idx_m,m_idx_m,dir_idx_m]=LME_var_linearTrack_3xtg_gen({k1,k2,k3,k4},all_infoscore_dir1,all_infoscore_dir2,pc1,pc2,1,AD_idx);

type_idx_m1=type_idx_m; % ntg old
type_idx_m1(type_idx_m1==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

%% sparsity LME
% suppose the k1-k4 as already been calculated

[var,trial_idx_m,type_idx_m,m_idx_m,dir_idx_m]=LME_var_linearTrack_3xtg_gen({k1,k2,k3,k4},spr1,spr2,[],[],1,AD_idx);

type_idx_m1=type_idx_m; % ntg old
type_idx_m1(type_idx_m1==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

type_idx_m1=type_idx_m;
type_idx_m1(type_idx_m1==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx_m1)},{categorical(trial_idx_m),categorical(m_idx_m),categorical(dir_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)+(fix1|random3)');

