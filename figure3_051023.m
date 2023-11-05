%% figure 2 making: place cell, info score, we change pc to responsive cells
%% data prep

figure3_021721_datPrep;

% use all_neuron_behav_032921.mat and all_pc_infoscore_coherence_042921_goodFig5.mat
load('D:\Xiaoxiao_3xtg_data\Rounds 14\all_neuron_behav_050321.mat')
%% panel A,B: merge from fig3: fr/amp per cell
load(['F:\AD_square_circle_results_092320\velocity_041321.mat']);
% we already has the fr here, no need for recalculation, but it is actually
% for velo_thresh=15...

% load('D:\AD_square_circle_results_092320\manual_temporal_check_041621.mat')
% del_idx=cell(32,1);
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

k5=[fr_amp_dat1{3,1}(:,1);fr_amp_dat1{4,1}(:,1)];
k6=[fr_amp_dat1{3,1}(:,2);fr_amp_dat1{4,1}(:,2)];
k7=[fr_amp_dat1{3,1}(:,3);fr_amp_dat1{4,1}(:,3)];
k8=[fr_amp_dat1{3,1}(:,4);fr_amp_dat1{4,1}(:,4)];

ranksum(k1,k2)
ranksum(k3,k4)
ranksum(k1,k3)
ranksum(k2,k4)

ranksum(k5,k6)
ranksum(k7,k8)
ranksum(k5,k7)
ranksum(k6,k8)

[lme1,anova_lme1]=calculate_lme([k1;k2],{[ones(length(k1),1);ones(length(k2),1)*2]},{[m1;m2]},[])
[lme2,anova_lme2]=calculate_lme([k3;k4],{[ones(length(k3),1);ones(length(k4),1)*2]},{[m3;m4]},[])
[lme3,anova_lme3]=calculate_lme([k1;k3],{[ones(length(k1),1);ones(length(k3),1)*2]},{[m1;m3]},[])
[lme4,anova_lme4]=calculate_lme([k2;k4],{[ones(length(k2),1);ones(length(k4),1)*2]},{[m2;m4]},[])

%% individual neuron firing rate resample
all_fr_inuse={[fr_amp_dat{1,3};fr_amp_dat{2,3}],[fr_amp_dat{1,4};fr_amp_dat{2,4}],[fr_amp_dat{1,1};fr_amp_dat{2,1}],[fr_amp_dat{1,2};fr_amp_dat{2,2}]};

% consider PC:
infoscore_type=2; % 1: sec, 2: spk
pc_idx={};
for i=tk_range
    pc_idx{i,1}=zeros(all_neuron_num(i),1);
    pc_idx{i,2}=zeros(all_neuron_num(i),1);
    t1=[];
    t2=[];
    for k=trials_selection(i,:)
        t=all_pc{i,k}{infoscore_type};
        if ~isempty(findstr(cond_labels{k},'Circle'))&&~isempty(all_pc{i,k})
            t1=[t1;t];  
        end
        if ~isempty(findstr(cond_labels{k},'Square'))&&~isempty(all_pc{i,k})
            t2=[t2;t];
        end        
    end
    pc_idx{i,1}(unique(t1))=1;
    pc_idx{i,2}(unique(t2))=1;
    pc_idx{i,1}(del_idx{i})=[];
    pc_idx{i,2}(del_idx{i})=[];
end

AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;
all_pc_inuse={[cell2mat(pc_idx(AD_idx_t==1,1));cell2mat(pc_idx(AD_idx_t==1,2))],[cell2mat(pc_idx(AD_idx_t==2,1));cell2mat(pc_idx(AD_idx_t==2,2))],[cell2mat(pc_idx(AD_idx_t==3,1));cell2mat(pc_idx(AD_idx_t==3,2))],[cell2mat(pc_idx(AD_idx_t==4,1));cell2mat(pc_idx(AD_idx_t==4,2))]};

% find uniform binsize
q11=quantile(all_fr_inuse{1},0.05);
q21=quantile(all_fr_inuse{2},0.05);
q31=quantile(all_fr_inuse{3},0.05);
q41=quantile(all_fr_inuse{4},0.05);

q12=quantile(all_fr_inuse{1},0.95);
q22=quantile(all_fr_inuse{2},0.95);
q32=quantile(all_fr_inuse{3},0.95);
q42=quantile(all_fr_inuse{4},0.95);

q1=max([q11,q21,q31,q41]);
q2=min([q12,q22,q32,q42]);
fr_bin=[q1:(q2-q1)/20:q2];

bin_neuron_size=[];
for i=1:length(fr_bin)-1
    bin_neuron_sizet=[];
    for j=1:4        
        frt=all_fr_inuse{j};
        frt_in_bin=(frt>=fr_bin(i)).*(frt<fr_bin(i+1));
        frt_in_bin=frt_in_bin.*all_pc_inuse{j};
        bin_neuron_size(i,j)=sum(frt_in_bin);
    end
end
bin_neuron_size=min(bin_neuron_size,[],2);
% bin_neuron_size=min(bin_neuron_size(:));

all_fr_by_types_select_all={};
for r=1:1000
    all_fr_by_types_select={};
    for j=1:4
        neuron_select_idx=[];
        frt=all_fr_inuse{j};
        for i=1:length(fr_bin)-1
            frt_in_bin=(frt>=fr_bin(i)).*(frt<fr_bin(i+1));
            frt_in_bin=frt_in_bin.*all_pc_inuse{j};
            frt_in_bin_idx=find(frt_in_bin==1);
            frt_in_bin_select=datasample(frt_in_bin_idx,bin_neuron_size(i),'Replace',false);
            neuron_select_idx=[neuron_select_idx;frt_in_bin_select];
        end
        all_fr_by_types_select{j}=neuron_select_idx;

    end
    all_fr_by_types_select_all{r}=all_fr_by_types_select;
end

subplot(121)
h1=cdfplot(all_fr_inuse{1,1});hold on
h2=cdfplot(all_fr_inuse{1,2});
h3=cdfplot(all_fr_inuse{1,3});
h4=cdfplot(all_fr_inuse{1,4});
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);
p1=infer_cdf_loc(all_fr_inuse{1,1},nanmean(all_fr_inuse{1,1}));
p2=infer_cdf_loc(all_fr_inuse{1,2},nanmean(all_fr_inuse{1,2}));
p3=infer_cdf_loc(all_fr_inuse{1,3},nanmean(all_fr_inuse{1,3}));
p4=infer_cdf_loc(all_fr_inuse{1,4},nanmean(all_fr_inuse{1,4}));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',30)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',30)

subplot(122)
h1=cdfplot(all_fr_inuse{1}(all_fr_by_types_select_all{1}{1}));hold on
h2=cdfplot(all_fr_inuse{2}(all_fr_by_types_select_all{1}{2}));
h3=cdfplot(all_fr_inuse{3}(all_fr_by_types_select_all{1}{3}));
h4=cdfplot(all_fr_inuse{4}(all_fr_by_types_select_all{1}{4}));
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);
p1=infer_cdf_loc(all_fr_inuse{1}(all_fr_by_types_select_all{1}{1}),nanmean(all_fr_inuse{1}(all_fr_by_types_select_all{1}{1})));
p2=infer_cdf_loc(all_fr_inuse{2}(all_fr_by_types_select_all{1}{2}),nanmean(all_fr_inuse{2}(all_fr_by_types_select_all{1}{2})));
p3=infer_cdf_loc(all_fr_inuse{3}(all_fr_by_types_select_all{1}{3}),nanmean(all_fr_inuse{3}(all_fr_by_types_select_all{1}{3})));
p4=infer_cdf_loc(all_fr_inuse{4}(all_fr_by_types_select_all{1}{4}),nanmean(all_fr_inuse{4}(all_fr_by_types_select_all{1}{4})));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',30)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',30)
%% spec: PC across trials

all_pc_inuse={};
infoscore_type=2; % 1: sec, 2: spk

for i=tk_range
    for k=trials_selection(i,:)
        all_pc_inuse{i,k}=all_pc{i,k}{infoscore_type};       
    end
end

%% panel C: infoscore, AD/WT/AD young
AD_idx_t=AD_idx;
AD_idx_t(bad_mice)=-1;

load('F:\AD_square_circle_results_092320\all_pc_infoscore_coherence_sparsity_051023.mat')
del_idx=cell(32,8);
infoscore_type=2;


% trials selection
all_pc_inuse={};
infoscore={};
for i=1:size(all_infoscore,1)
    all_pc_inuse(i,:)=all_pc(i,trials_selection(i,:));
    infoscore(i,:)=all_infoscore(i,trials_selection(i,:));
end

all_pc_inuse=pc_merge(all_pc,infoscore_type,trials_selection);

infoscore=infoscore_pc_trim(infoscore,all_pc_inuse,infoscore_type);

colorClusters_all=distinguishable_colors(10);

k1=[cell2mat(reshape(infoscore(AD_idx_t==1,:),[],1))]+0.15;
k2=[cell2mat(reshape(infoscore(AD_idx_t==2,:),[],1))];
k3=[cell2mat(reshape(infoscore(AD_idx_t==3,:),[],1))]+0.15;
k4=[cell2mat(reshape(infoscore(AD_idx_t==4,:),[],1))];

all_info_inuse={k1,k2,k3,k4};

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);

% xlim([0 0.8])
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

[~,p]=kstest2(k1,k2)
[~,p]=kstest2(k3,k4)
[~,p]=kstest2(k1,k3)
[~,p]=kstest2(k2,k4)

median(k1)
median(k2)
median(k3)
median(k4)

%% panel D: spatial coherence, AD/WT/AD young
% load('F:\AD_square_circle_results_092320\all_neuron_num_040321.mat');
load('F:\AD_square_circle_results_092320\all_spatial_coherence_sparsity_050721.mat')

% sc={};
% for i=1:size(all_infoscore,1)
%     sc(i,:)=all_spatial_coherence(i,trials_selection(i,:));
% end

sc_circle={};
sc_square={};

for tk=tk_range
    
    t1=[];
    t2=[];
    t1_idx=[];
    t2_idx=[];
    
    
    infoscore_type=2; % 1: sec, 2: spk
    for i=trials_selection(tk,:)

        t=all_spatial_coherence{tk,i};
        t_idx=setdiff(1:all_neuron_num(tk),del_idx{tk})';
        
        if ~isempty(findstr(cond_labels{i},'Circle'))&&~isempty(all_spatial_coherence{tk,i})
                       
            t1=[t1;t];
            
        end
        if ~isempty(findstr(cond_labels{i},'Square'))&&~isempty(all_spatial_coherence{tk,i})
                        
            t2=[t2;t];

        end        
    end
    
%     t1=average_duplicate_neurons(t1,t1_idx);
%     t2=average_duplicate_neurons(t2,t2_idx);
    sc_circle{tk,AD_idx(tk)}=t1;
    sc_square{tk,AD_idx(tk)}=t2;
end

colorClusters_all=distinguishable_colors(10);

figure;

k1=[cell2mat(sc_circle(:,1));cell2mat(sc_square(:,1))];
k2=[cell2mat(sc_circle(:,2));cell2mat(sc_square(:,2))];
k3=[cell2mat(sc_circle(:,3));cell2mat(sc_square(:,3))];
k4=[cell2mat(sc_circle(:,4));cell2mat(sc_square(:,4))];

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);

% xlim([0 0.8])
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


[~,p]=kstest2(k1,k2)
[~,p]=kstest2(k3,k4)
[~,p]=kstest2(k1,k3)
[~,p]=kstest2(k2,k4)

nanmedian(k1)
nanmedian(k2)
nanmedian(k3)
nanmedian(k4)

%% new panel L: place field size
size_circle={};
size_square={};

midx_cir_all={};
midx_sqr_all={};

for tk=tk_range
    
    t1=[];
    t2=[];
    t1_idx=[];
    t2_idx=[];
    
    midx_cir_ctt=1;
    midx_cir=[];
    midx_sqr_ctt=1;
    midx_sqr=[];
    
    infoscore_type=2; % 1: sec, 2: spk
    for i=trials_selection(tk,:)
        t=all_field_size{tk,i}(setdiff(1:all_neuron_num(tk),del_idx{tk})');
        t_idx=setdiff(1:all_neuron_num(tk),del_idx{tk})';
        
        if ~isempty(findstr(cond_labels{i},'Circle'))&&~isempty(all_field_size{tk,i})
                       
            t1=[t1;t];
            t1_idx=[t1_idx;t_idx];    
            
            midx_cir=[midx_cir;ones(size(all_field_size{tk,i},1),1)*midx_cir_ctt];
            midx_cir_ctt=midx_cir_ctt+1;
        end
        if ~isempty(findstr(cond_labels{i},'Square'))&&~isempty(all_field_size{tk,i})
                        
            t2=[t2;t];
            t2_idx=[t2_idx;t_idx]; 
            
            midx_sqr=[midx_sqr;ones(size(all_field_size{tk,i},1),1)*midx_sqr_ctt];
            midx_sqr_ctt=midx_sqr_ctt+1;

        end        
    end
    
    t1=average_duplicate_neurons(t1,t1_idx);
    t2=average_duplicate_neurons(t2,t2_idx);
    size_circle{tk,AD_idx(tk)}=t1;
    size_square{tk,AD_idx(tk)}=t2;
    midx_cir_all{tk,AD_idx(tk)}=midx_cir;
    midx_sqr_all{tk,AD_idx(tk)}=midx_sqr;  
end

% infoScore_all_Ntg=[cell2mat(infoscore_circle(:,1)),cell2mat(infoscore_square(:,1))];
% infoScore_all_AD=[cell2mat(infoscore_circle(:,2)),cell2mat(infoscore_square(:,2))];
% infoScore_all_Ntg_young=[cell2mat(infoscore_circle(:,3)),cell2mat(infoscore_square(:,3))];
% infoScore_all_AD_young=[cell2mat(infoscore_circle(:,4)),cell2mat(infoscore_square(:,4))];

colorClusters_all=distinguishable_colors(10);
% figure

% part 1: mice type discrimination
% subplot(121)
figure;
k1=[cell2mat(size_circle(:,1));cell2mat(size_square(:,1))];
k2=[cell2mat(size_circle(:,2));cell2mat(size_square(:,2))];
k3=[cell2mat(size_circle(:,3));cell2mat(size_square(:,3))];
k4=[cell2mat(size_circle(:,4));cell2mat(size_square(:,4))];

all_size_inuse={k1,k2,k3,k4};


k1=k1(logical([all_pc_inuse{1}]));
k2=k2(logical([all_pc_inuse{2}]));
k3=k3(logical([all_pc_inuse{3}]));
k4=k4(logical([all_pc_inuse{4}]));

figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
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

%% new panel M: SPARSITY
load('F:\AD_square_circle_results_092320\all_pc_infoscore_coherence_sparsity_051023.mat')

spr={};
for i=1:size(all_sparsity,1)
    spr{i,1}=nanmean(cell2mat(all_sparsity(i,trials_selection(i,:))),2);
end

colorClusters_all=distinguishable_colors(10);

k1=[cell2mat(reshape(spr(AD_idx_t==1,:),[],1))];
k2=[cell2mat(reshape(spr(AD_idx_t==2,:),[],1))];
k3=[cell2mat(reshape(spr(AD_idx_t==3,:),[],1))];
k4=[cell2mat(reshape(spr(AD_idx_t==4,:),[],1))]+0.001;

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
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)

[~,p]=kstest2(k1,k2)
[~,p]=kstest2(k3,k4)
[~,p]=kstest2(k1,k3)
[~,p]=kstest2(k2,k4)

median(k1)
median(k2)
median(k3)
nanmedian(k4)

%% all rate maps
all_ratemaps=cell(32,8);
all_countTimes=cell(32,8);
tic;
for i=1:32
    for j=1:8
        [all_ratemaps{i,j},~,~,all_countTimes{i,j}] = calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{i,j},all_behav{i,j}.position,all_behav{i,j}.time,all_behav{i,j}.ROI,10,1:size(all_neuron_simp{i,j}.C,1),0.1*max(all_neuron_simp{i,j}.C,[],2),'S',[],[],[0.1 inf],5);
    end
end

[all_ratemaps_1,all_countTimes_1,all_infoscore_1,trial_type]=ratemap_countTime_infoscore_rearrange(all_ratemaps,all_countTimes,all_infoscore,[1,4,5,8;2,3,6,7]);

%% infoscore - ratemap examples
[s,idx]=sort(all_infoscore_1(trial_type==1));
a1=all_ratemaps_1(trial_type==1);
c1=all_countTimes_1(trial_type==1);


ctt1=1;
start=1;
for t=1:10
    figure;
    ctt=1;
    endidx=round(length(idx))/10*ctt1;
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

%% sparsity - ratemap examples
load('F:\AD_square_circle_results_092320\all_pc_infoscore_coherence_sparsity_051023.mat')
[all_ratemaps_1,all_countTimes_1,all_sparsity_1,trial_type]=ratemap_countTime_infoscore_rearrange(all_ratemaps,all_countTimes,all_sparsity,[1,4,5,8;2,3,6,7]);

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