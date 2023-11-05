% figure 4 panel A
AD_idx_t=AD_idx;

load(['D:\AD_square_circle_results_092320\velocity_041621.mat']);
load('D:\AD_square_circle_results_092320\manual_temporal_check_041621.mat')
% we already has the fr here, no need for recalculation, but it is actually
% for velo_thresh=15...
% 3: trim neuronIndividuals_new
fr_cir_run={};
fr_cir_nrun={};
fr_sqr_run={};
fr_sqr_nrun={};

amp_cir_run={};
amp_cir_nrun={};
amp_sqr_run={};
amp_sqr_nrun={};

fr_cir_run_chk={};
fr_cir_nrun_chk={};

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
    
    [fr_cir_run{tk},amp_cir_run{tk},fr_cir_run_chk{tk}]=general_fr_amp_calculation(all_neuron_run(tk,trial_select_cir),del_idx{tk},responsive_thresh); % neuron trim done here 0.02Hz can be applied locally
    [fr_cir_nrun{tk},amp_cir_nrun{tk},fr_cir_nrun_chk{tk}]=general_fr_amp_calculation(all_neuron_nrun(tk,trial_select_cir),del_idx{tk},responsive_thresh);
    [fr_sqr_run{tk},amp_sqr_run{tk}]=general_fr_amp_calculation(all_neuron_run(tk,trial_select_sqr),del_idx{tk},responsive_thresh);
    [fr_sqr_nrun{tk},amp_sqr_nrun{tk}]=general_fr_amp_calculation(all_neuron_nrun(tk,trial_select_sqr),del_idx{tk},responsive_thresh);

end
fr_amp_dat={};
fr_amp_dat{1,1}=cell2mat(fr_cir_run(:,AD_idx==3)');
fr_amp_dat{1,2}=cell2mat(fr_cir_run(:,AD_idx==4)');
fr_amp_dat{1,3}=cell2mat(fr_cir_run(:,AD_idx==1)');
fr_amp_dat{1,4}=cell2mat(fr_cir_run(:,AD_idx==2)');

fr_amp_dat{1,5}=cell2mat(fr_cir_nrun(:,AD_idx==3)');
fr_amp_dat{1,6}=cell2mat(fr_cir_nrun(:,AD_idx==4)');
fr_amp_dat{1,7}=cell2mat(fr_cir_nrun(:,AD_idx==1)');
fr_amp_dat{1,8}=cell2mat(fr_cir_nrun(:,AD_idx==2)');

fr_amp_dat{2,1}=cell2mat(fr_sqr_run(:,AD_idx==3)');
fr_amp_dat{2,2}=cell2mat(fr_sqr_run(:,AD_idx==4)');
fr_amp_dat{2,3}=cell2mat(fr_sqr_run(:,AD_idx==1)');
fr_amp_dat{2,4}=cell2mat(fr_sqr_run(:,AD_idx==2)');

fr_amp_dat{2,5}=cell2mat(fr_sqr_nrun(:,AD_idx==3)');
fr_amp_dat{2,6}=cell2mat(fr_sqr_nrun(:,AD_idx==4)');
fr_amp_dat{2,7}=cell2mat(fr_sqr_nrun(:,AD_idx==1)');
fr_amp_dat{2,8}=cell2mat(fr_sqr_nrun(:,AD_idx==2)');

fr_amp_dat{3,1}=cell2mat(amp_cir_run(:,AD_idx==3)');
fr_amp_dat{3,2}=cell2mat(amp_cir_run(:,AD_idx==4)');
fr_amp_dat{3,3}=cell2mat(amp_cir_run(:,AD_idx==1)');
fr_amp_dat{3,4}=cell2mat(amp_cir_run(:,AD_idx==2)');

fr_amp_dat{3,5}=cell2mat(amp_cir_nrun(:,AD_idx==3)');
fr_amp_dat{3,6}=cell2mat(amp_cir_nrun(:,AD_idx==4)');
fr_amp_dat{3,7}=cell2mat(amp_cir_nrun(:,AD_idx==1)');
fr_amp_dat{3,8}=cell2mat(amp_cir_nrun(:,AD_idx==2)');

fr_amp_dat{4,1}=cell2mat(amp_sqr_run(:,AD_idx==3)');
fr_amp_dat{4,2}=cell2mat(amp_sqr_run(:,AD_idx==4)');
fr_amp_dat{4,3}=cell2mat(amp_sqr_run(:,AD_idx==1)');
fr_amp_dat{4,4}=cell2mat(amp_sqr_run(:,AD_idx==2)');

fr_amp_dat{4,5}=cell2mat(amp_sqr_nrun(:,AD_idx==3)');
fr_amp_dat{4,6}=cell2mat(amp_sqr_nrun(:,AD_idx==4)');
fr_amp_dat{4,7}=cell2mat(amp_sqr_nrun(:,AD_idx==1)');
fr_amp_dat{4,8}=cell2mat(amp_sqr_nrun(:,AD_idx==2)');


for i=1:2
    for j=1:size(fr_amp_dat,2)
        fr_amp_dat{i,j}(fr_amp_dat{i,j}<=responsive_thresh)=[];
        fr_amp_dat{i,j}(fr_amp_dat{i,j}>1)=[];
    end
end
for i=3:4
    for j=1:size(fr_amp_dat,2)
        fr_amp_dat{i-2,j}(fr_amp_dat{i-2,j}<=responsive_thresh)=[];
        fr_amp_dat{i,j}(fr_amp_dat{i,j}>150)=[];
    end
end

% 
fr_amp_dat1=fill_nan_to_cellmat(fr_amp_dat) %T first four column: run

k1=[fr_amp_dat1{1,1}(:,1);fr_amp_dat1{2,1}(:,1)];
k2=[fr_amp_dat1{1,1}(:,2);fr_amp_dat1{2,1}(:,2)];
k3=[fr_amp_dat1{1,1}(:,3);fr_amp_dat1{2,1}(:,3)];
k4=[fr_amp_dat1{1,1}(:,4);fr_amp_dat1{2,1}(:,4)];


type_cir_idx={};
trial_cir_idx={};
m_cir_idx={};

type_sqr_idx={};
trial_sqr_idx={};
m_sqr_idx={};

for i=tk_range
    type_cir_idx{i,1}=ones(size(fr_cir_run{i},1),1)*AD_idx(i);
    trial_cir_idx{i,1}=ones(size(fr_cir_run{i},1),1);
    m_cir_idx{i,1}=ones(size(fr_cir_run{i},1),1)*i;
    
    type_sqr_idx{i,1}=ones(size(fr_cir_run{i},1),1)*AD_idx(i);
    trial_sqr_idx{i,1}=ones(size(fr_cir_run{i},1),1);
    m_sqr_idx{i,1}=ones(size(fr_cir_run{i},1),1)*i;
end

trial_cir_idx_all={cell2mat(type_cir_idx(AD_idx==1)),cell2mat(type_cir_idx(AD_idx==2)),cell2mat(type_cir_idx(AD_idx==3)),cell2mat(type_cir_idx(AD_idx==4))};
trial_sqr_idx_all={cell2mat(type_sqr_idx(AD_idx==1)),cell2mat(type_sqr_idx(AD_idx==2)),cell2mat(type_sqr_idx(AD_idx==3)),cell2mat(type_sqr_idx(AD_idx==4))};

type_cir_idx_all={cell2mat(trial_cir_idx(AD_idx==1)),cell2mat(trial_cir_idx(AD_idx==2)),cell2mat(trial_cir_idx(AD_idx==3)),cell2mat(trial_cir_idx(AD_idx==4))};
type_sqr_idx_all={cell2mat(trial_sqr_idx(AD_idx==1)),cell2mat(trial_sqr_idx(AD_idx==2)),cell2mat(trial_sqr_idx(AD_idx==3)),cell2mat(trial_sqr_idx(AD_idx==4))};

m_cir_idx_all={cell2mat(m_cir_idx(AD_idx==1)),cell2mat(m_cir_idx(AD_idx==2)),cell2mat(m_cir_idx(AD_idx==3)),cell2mat(m_cir_idx(AD_idx==4))};
m_sqr_idx_all={cell2mat(m_sqr_idx(AD_idx==1)),cell2mat(m_sqr_idx(AD_idx==2)),cell2mat(m_sqr_idx(AD_idx==3)),cell2mat(m_sqr_idx(AD_idx==4))};

trial_cir_idx_all_m=fill_nan_to_cellmat(trial_cir_idx_all);
trial_sqr_idx_all_m=fill_nan_to_cellmat(trial_sqr_idx_all);
type_cir_idx_all_m=fill_nan_to_cellmat(type_cir_idx_all);
type_sqr_idx_all_m=fill_nan_to_cellmat(type_sqr_idx_all);
m_cir_idx_all_m=fill_nan_to_cellmat(m_cir_idx_all);
m_sqr_idx_all_m=fill_nan_to_cellmat(m_sqr_idx_all);

var=[k1;k2;k3;k4];
% trial_idx_m=[[trial_cir_idx_all_m{1}(:,1);trial_sqr_idx_all_m{1}(:,1)];[trial_cir_idx_all_m{1}(:,2);trial_sqr_idx_all_m{1}(:,2)];[trial_cir_idx_all_m{1}(:,3);trial_sqr_idx_all_m{1}(:,3)];[trial_cir_idx_all_m{1}(:,4);trial_sqr_idx_all_m{1}(:,4)]];
% type_idx_m=[[type_cir_idx_all_m{1}(:,1);type_sqr_idx_all_m{1}(:,1)];[type_cir_idx_all_m{1}(:,2);type_sqr_idx_all_m{1}(:,2)];[type_cir_idx_all_m{1}(:,3);type_sqr_idx_all_m{1}(:,3)];[type_cir_idx_all_m{1}(:,4);type_sqr_idx_all_m{1}(:,4)]];
% m_idx_m=[[m_cir_idx_all_m{1}(:,1);m_sqr_idx_all_m{1}(:,1)];[m_cir_idx_all_m{1}(:,2);m_sqr_idx_all_m{1}(:,2)];[m_cir_idx_all_m{1}(:,3);m_sqr_idx_all_m{1}(:,3)];[m_cir_idx_all_m{1}(:,4);m_sqr_idx_all_m{1}(:,4)]];
trial_idx_m=[ones(length(k1)/2,1);ones(length(k1)/2,1)*2;ones(length(k2)/2,1);ones(length(k2)/2,1)*2;ones(length(k3)/2,1);ones(length(k3)/2,1)*2;ones(length(k4)/2,1);ones(length(k4)/2,1)*2];
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[[m_cir_idx_all_m{1}(:,1);m_sqr_idx_all_m{1}(:,1)];[m_cir_idx_all_m{1}(:,2);m_sqr_idx_all_m{1}(:,2)];[m_cir_idx_all_m{1}(:,3);m_sqr_idx_all_m{1}(:,3)];[m_cir_idx_all_m{1}(:,4);m_sqr_idx_all_m{1}(:,4)]];

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

ranksum(k1,k2) % young ntg VS YOUNG ad, although for lme the ad_idx is corrected, it is not the case for k1-k4 
ranksum(k3,k4) % old ntg vs old AD 
ranksum(k1,k3)
ranksum(k2,k4)
%% panel B amplitude

k5=[fr_amp_dat1{3,1}(:,1);fr_amp_dat1{4,1}(:,1)];
k6=[fr_amp_dat1{3,1}(:,2);fr_amp_dat1{4,1}(:,2)];
k7=[fr_amp_dat1{3,1}(:,3);fr_amp_dat1{4,1}(:,3)];
k8=[fr_amp_dat1{3,1}(:,4);fr_amp_dat1{4,1}(:,4)];

var=[k5;k6;k7;k8];
% trial_idx_m=[[trial_cir_idx_all_m{1}(:,1);trial_sqr_idx_all_m{1}(:,1)];[trial_cir_idx_all_m{1}(:,2);trial_sqr_idx_all_m{1}(:,2)];[trial_cir_idx_all_m{1}(:,3);trial_sqr_idx_all_m{1}(:,3)];[trial_cir_idx_all_m{1}(:,4);trial_sqr_idx_all_m{1}(:,4)]];
% type_idx_m=[[type_cir_idx_all_m{1}(:,1);type_sqr_idx_all_m{1}(:,1)];[type_cir_idx_all_m{1}(:,2);type_sqr_idx_all_m{1}(:,2)];[type_cir_idx_all_m{1}(:,3);type_sqr_idx_all_m{1}(:,3)];[type_cir_idx_all_m{1}(:,4);type_sqr_idx_all_m{1}(:,4)]];
% m_idx_m=[[m_cir_idx_all_m{1}(:,1);m_sqr_idx_all_m{1}(:,1)];[m_cir_idx_all_m{1}(:,2);m_sqr_idx_all_m{1}(:,2)];[m_cir_idx_all_m{1}(:,3);m_sqr_idx_all_m{1}(:,3)];[m_cir_idx_all_m{1}(:,4);m_sqr_idx_all_m{1}(:,4)]];
trial_idx_m=[ones(length(k1)/2,1);ones(length(k1)/2,1)*2;ones(length(k2)/2,1);ones(length(k2)/2,1)*2;ones(length(k3)/2,1);ones(length(k3)/2,1)*2;ones(length(k4)/2,1);ones(length(k4)/2,1)*2];
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[[m_cir_idx_all_m{1}(:,1);m_sqr_idx_all_m{1}(:,1)];[m_cir_idx_all_m{1}(:,2);m_sqr_idx_all_m{1}(:,2)];[m_cir_idx_all_m{1}(:,3);m_sqr_idx_all_m{1}(:,3)];[m_cir_idx_all_m{1}(:,4);m_sqr_idx_all_m{1}(:,4)]];

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

ranksum(k5,k6)
ranksum(k7,k8)
ranksum(k5,k7)
ranksum(k6,k8)
%% decoding LME (per trial)
AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;

type_idx={};
trial_idx={};
m_idx={};
for i=tk_range
    type_idx{i,1}=[AD_idx(i),AD_idx(i),AD_idx(i),AD_idx(i)]
    trial_idx{i,1}=trials_selection(i,:);
    m_idx{i,1}=[i,i,i,i];
end

decoded_accuracy_allTrial_cell=[reshape(cell2mat(decoded_accuracy(AD_idx_t==1)),[],1);reshape(cell2mat(decoded_accuracy(AD_idx_t==2)),[],1);reshape(cell2mat(decoded_accuracy(AD_idx_t==3)),[],1);reshape(cell2mat(decoded_accuracy(AD_idx_t==4)),[],1)];
trial_idx_m=[reshape(cell2mat(trial_idx(AD_idx_t==1)),[],1);reshape(cell2mat(trial_idx(AD_idx_t==2)),[],1);reshape(cell2mat(trial_idx(AD_idx_t==3)),[],1);reshape(cell2mat(trial_idx(AD_idx_t==4)),[],1)];
type_idx_m=[reshape(cell2mat(type_idx(AD_idx_t==1)),[],1);reshape(cell2mat(type_idx(AD_idx_t==2)),[],1);reshape(cell2mat(type_idx(AD_idx_t==3)),[],1);reshape(cell2mat(type_idx(AD_idx_t==4)),[],1)];
m_idx_m=[reshape(cell2mat(m_idx(AD_idx_t==1)),[],1);reshape(cell2mat(m_idx(AD_idx_t==2)),[],1);reshape(cell2mat(m_idx(AD_idx_t==3)),[],1);reshape(cell2mat(m_idx(AD_idx_t==4)),[],1)];

type_decode=type_idx_m;
type_decode(type_decode==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(decoded_accuracy_allTrial_cell,{categorical(type_decode)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_decode=type_idx_m;
type_decode(type_decode==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(decoded_accuracy_allTrial_cell,{categorical(type_decode)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_decode=type_idx_m;
type_decode(type_decode==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(decoded_accuracy_allTrial_cell,{categorical(type_decode)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_decode=type_idx_m;
type_decode(type_decode==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(decoded_accuracy_allTrial_cell,{categorical(type_decode)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');



%% infoscore
AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;
[var,trial_idx_m,type_idx_m,m_idx_m]=LME_var_openField_3xtg_gen({k1,k2,k3,k4},all_infoscore,all_pc_inuse,trials_selection,AD_idx_t);

type_idx=type_idx_m;
type_idx(type_idx==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

%% coherence
AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;
[var,trial_idx_m,type_idx_m,m_idx_m]=LME_var_openField_3xtg_gen({k1,k2,k3,k4},all_spatial_coherence,[],trials_selection,AD_idx_t);

type_idx=type_idx_m;
type_idx(type_idx==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');


%% sparsity
AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;
[var,trial_idx_m,type_idx_m,m_idx_m]=LME_var_openField_3xtg_gen({k1,k2,k3,k4},spr,[],trials_selection,AD_idx_t);

type_idx=type_idx_m;
type_idx(type_idx==1)=0; % reference group is 1 or 0
[lme1,anova_lme1]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==2)=0; % reference group is 1 or 0
[lme2,anova_lme2]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==3)=0; % reference group is 1 or 0
[lme3,anova_lme3]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');

type_idx=type_idx_m;
type_idx(type_idx==4)=0; % reference group is 1 or 0
[lme4,anova_lme4]=calculate_lme_011221(var,{categorical(type_idx)},{categorical(trial_idx_m),categorical(m_idx_m)},'y~fix1+(fix1|random1)+(fix1|random2)');
