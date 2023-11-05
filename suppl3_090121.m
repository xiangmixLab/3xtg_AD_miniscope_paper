%% suppl 3: activity coverage

%% OF activate bin fractions
load('F:\AD_square_circle_results_092320\all_neuron_behav_032921.mat');
responsive_thresh=0;
del_idx=cell(32,1);
load('F:\AD_square_circle_results_092320\exp_info_etgAD_031621.mat');
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


tic;
binsize=10;
fr_all={};
count_all={};
ctime={};
for tk=1:size(all_neuron_simp,1)
%     load([folderName{tk},'\','neuronIndividuals_new.mat']);
    for j=1:size(all_neuron_simp,2)
        [fr_all{tk,j},count_all{tk,j},~,ctime{tk,j}] = calculatingCellSpatialForSingleData_Suoqin(all_neuron_simp{tk,j},all_behav{tk,j}.position,all_behav{tk,j}.time,[0 0 max(all_behav{tk,j}.position(:,1)),max(all_behav{tk,j}.position(:,2))],binsize,1:size(all_neuron_simp{tk,j}.C,1),3*std(all_neuron_simp{tk,j}.S,[],2),'S',[],[],[0.1 inf],5);

    end
    toc;
end


AD_idx=[4 4 4 4 3 3 3 3 4 4 3 3 4 4 2 2 2 1 1 1 2 2 1 2 1 1 1 2 2 1 1 1];
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

fraction_all={};
midx={};
for tk=tk_range
    fraction_all{tk}=[];
    midx{tk}=[];
    
    trial_select_cir=[];
    trial_select_sqr=[];
    for j=1:length(trials_selection(tk,:))
        if ~isempty(strfind(cond_labels{trials_selection(tk,j)},'Circle'))
            trial_select_cir=[trial_select_cir,trials_selection(tk,j)];
        else
            trial_select_sqr=[trial_select_sqr,trials_selection(tk,j)];
        end
    end
    
    midx_temp=[];
    fraction_cell_temp={};
    for j=trials_selection(tk,:)
        val=[];
        
        for k=1:length(fr_all{tk,j})            
            val(k,1) = nansum(nansum(fr_all{tk,j}{k}>0))/nansum(nansum(ctime{tk,j}>0));
            midx_temp(k,1)=tk;
        end
        fraction_cell_temp{tk,j}=[val];
    end
    fraction_all{tk}=[mean([fraction_cell_temp{tk,trial_select_cir(1)},fraction_cell_temp{tk,trial_select_cir(2)}],2);mean([fraction_cell_temp{tk,trial_select_sqr(1)},fraction_cell_temp{tk,trial_select_sqr(2)}],2)];
    midx{tk}=[midx_temp;midx_temp];
end

fraction_all_dat={};
fraction_all_dat{1,1}=cell2mat(fraction_all(:,AD_idx==3)');
fraction_all_dat{1,2}=cell2mat(fraction_all(:,AD_idx==4)');
fraction_all_dat{1,3}=cell2mat(fraction_all(:,AD_idx==1)');
fraction_all_dat{1,4}=cell2mat(fraction_all(:,AD_idx==2)');

fraction_all_dat1=fill_nan_to_cellmat(fraction_all_dat); %T first four column: run

% LME
k1=fraction_all_dat{1,1};
k2=fraction_all_dat{1,2};
k3=fraction_all_dat{1,3};
k4=fraction_all_dat{1,4};

m1=cell2mat(midx(:,AD_idx==3)');
m2=cell2mat(midx(:,AD_idx==4)');
m3=cell2mat(midx(:,AD_idx==1)');
m4=cell2mat(midx(:,AD_idx==2)');

var=[k1;k2;k3;k4];
trial_idx_m=[];
for i=1:4
    trial_idx_m=[trial_idx_m;ones(length(fraction_all_dat{1,i})/2,1);ones(length(fraction_all_dat{1,i})/2,1)*2];
end
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[m1;m2;m3;m4];

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

[mean(k1),mean(k2),mean(k3),mean(k4)]

%% OF bin rate deviation

rate_dev_all={};
midx={};
for tk=tk_range
    rate_dev_all{tk}=[];
    midx{tk}=[];
    
    trial_select_cir=[];
    trial_select_sqr=[];
    for j=1:length(trials_selection(tk,:))
        if ~isempty(strfind(cond_labels{trials_selection(tk,j)},'Circle'))
            trial_select_cir=[trial_select_cir,trials_selection(tk,j)];
        else
            trial_select_sqr=[trial_select_sqr,trials_selection(tk,j)];
        end
    end
    
    midx_temp=[];
    dev_cell_temp={};
    for j=trials_selection(tk,:)
        val=[];
        
        for k=1:length(fr_all{tk,j})            
            val(k,1) = nanmax(zscore(fr_all{tk,j}{k}(:)));
            midx_temp(k,1)=tk;
        end
        dev_cell_temp{tk,j}=[val];
    end
    rate_dev_all{tk}=[mean([dev_cell_temp{tk,trial_select_cir(1)},dev_cell_temp{tk,trial_select_cir(2)}],2);mean([dev_cell_temp{tk,trial_select_sqr(1)},dev_cell_temp{tk,trial_select_sqr(2)}],2)];
    midx{tk}=[midx_temp;midx_temp];
end

dev_all_dat={};
dev_all_dat{1,1}=cell2mat(rate_dev_all(:,AD_idx==3)');
dev_all_dat{1,2}=cell2mat(rate_dev_all(:,AD_idx==4)');
dev_all_dat{1,3}=cell2mat(rate_dev_all(:,AD_idx==1)');
dev_all_dat{1,4}=cell2mat(rate_dev_all(:,AD_idx==2)');

dev_all_dat1=fill_nan_to_cellmat(dev_all_dat); %T first four column: run

% LME
k1=fraction_all_dat{1,1};
k2=fraction_all_dat{1,2};
k3=fraction_all_dat{1,3};
k4=fraction_all_dat{1,4};

m1=cell2mat(midx(:,AD_idx==3)');
m2=cell2mat(midx(:,AD_idx==4)');
m3=cell2mat(midx(:,AD_idx==1)');
m4=cell2mat(midx(:,AD_idx==2)');

var=[k1;k2;k3;k4];
trial_idx_m=[];
for i=1:4
    trial_idx_m=[trial_idx_m;ones(length(fraction_all_dat{1,i})/2,1);ones(length(fraction_all_dat{1,i})/2,1)*2];
end
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[m1;m2;m3;m4];

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

[mean(k1),mean(k2),mean(k3),mean(k4)]

%% OF autocorr

autocorr_all={};
radius_all={};
localmax_all={};
field_example={};
midx={};
for tk=tk_range
    rate_dev_all{tk}=[];
    midx{tk}=[];
    
    trial_select_cir=[];
    trial_select_sqr=[];
    for j=1:length(trials_selection(tk,:))
        if ~isempty(strfind(cond_labels{trials_selection(tk,j)},'Circle'))
            trial_select_cir=[trial_select_cir,trials_selection(tk,j)];
        else
            trial_select_sqr=[trial_select_sqr,trials_selection(tk,j)];
        end
    end
    
    midx_temp=[];
    acorr_cell_temp={};
    radius_cell_temp={};
    localmax_cell_temp={};
    field_temp={};
    ctt=1;
    for j=trials_selection(tk,:)
        val=[];
        val2=[];
        val3=[];
        
        for k=1:length(fr_all{tk,j})            
            [~,~,~,val(k,1),val3(k,1),field_temp{k,ctt}]=firing_field_radius_autocorr(fr_all{tk,j}{k},1);
            [val2(k,1)]=firing_field_radius_autocorr(fr_all{tk,j}{k},1);
            midx_temp(k,1)=tk;
        end
        acorr_cell_temp{tk,j}=[val];
        radius_cell_temp{tk,j}=[val2];
        localmax_cell_temp{tk,j}=[val3];
        ctt=ctt+1;
    end
    autocorr_all{tk}=[nanmean([acorr_cell_temp{tk,trial_select_cir(1)},acorr_cell_temp{tk,trial_select_cir(2)}],2);nanmean([acorr_cell_temp{tk,trial_select_sqr(1)},acorr_cell_temp{tk,trial_select_sqr(2)}],2)];
    radius_all{tk}=[nanmean([radius_cell_temp{tk,trial_select_cir(1)},radius_cell_temp{tk,trial_select_cir(2)}],2);nanmean([radius_cell_temp{tk,trial_select_sqr(1)},radius_cell_temp{tk,trial_select_sqr(2)}],2)];
    localmax_all{tk}=[nanmean([localmax_cell_temp{tk,trial_select_cir(1)},localmax_cell_temp{tk,trial_select_cir(2)}],2);nanmean([localmax_cell_temp{tk,trial_select_sqr(1)},localmax_cell_temp{tk,trial_select_sqr(2)}],2)];
    
    field_example{tk}=reshape(field_temp,length(fr_all{tk,j})*length(trials_selection(tk,:)),1);
    midx{tk}=[midx_temp;midx_temp];
end


% following analysis
acorr_all_dat={};
acorr_all_dat{1,1}=cell2mat(autocorr_all(:,AD_idx==1)'); % Ntg old
acorr_all_dat{1,2}=cell2mat(autocorr_all(:,AD_idx==2)');
acorr_all_dat{1,3}=cell2mat(autocorr_all(:,AD_idx==3)');
acorr_all_dat{1,4}=cell2mat(autocorr_all(:,AD_idx==4)');

acorr_all_dat1=fill_nan_to_cellmat(acorr_all_dat([3,4,1,2])); %T first four column: run

radius_all_dat={};
radius_all_dat{1,1}=cell2mat(radius_all(:,AD_idx==1)'); % Ntg old
radius_all_dat{1,2}=cell2mat(radius_all(:,AD_idx==2)');
radius_all_dat{1,3}=cell2mat(radius_all(:,AD_idx==3)');
radius_all_dat{1,4}=cell2mat(radius_all(:,AD_idx==4)');

radius_all_dat1=fill_nan_to_cellmat(radius_all_dat([3,4,1,2])); %T first four column: run

localmax_all_dat={};
localmax_all_dat{1,1}=cell2mat(localmax_all(:,AD_idx==1)'); % Ntg old
localmax_all_dat{1,2}=cell2mat(localmax_all(:,AD_idx==2)');
localmax_all_dat{1,3}=cell2mat(localmax_all(:,AD_idx==3)');
localmax_all_dat{1,4}=cell2mat(localmax_all(:,AD_idx==4)');

localmax_all_dat1=fill_nan_to_cellmat(localmax_all_dat([3,4,1,2])); %T first four column: run

% area cdf
colorClusters_all=distinguishable_colors(10);
figure;
k1=acorr_all_dat{1,1};
k2=acorr_all_dat{1,2};
idx=find(k2>1000)
k2(idx)=[]; % artifact, only 1 val >1000 while others <300
k3=acorr_all_dat{1,3};
k4=acorr_all_dat{1,4};

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[ 0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);
p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',30)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',30)

[h,p1]=kstest2(k1,k2)
[h,p2]=kstest2(k3,k4)

% radius cdf
figure;
k1=radius_all_dat{1,1};
k2=radius_all_dat{1,2};
k2(idx)=[]; % artifact,
k3=radius_all_dat{1,3};
k4=radius_all_dat{1,4};

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[ 0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);
p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',30)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',30)

% localmax lme

k1=localmax_all_dat{1,3};
k2=localmax_all_dat{1,4};
k3=localmax_all_dat{1,1};
k4=localmax_all_dat{1,2};

m1=cell2mat(midx(:,AD_idx==3)');
m2=cell2mat(midx(:,AD_idx==4)');
m3=cell2mat(midx(:,AD_idx==1)');
m4=cell2mat(midx(:,AD_idx==2)');

var=[k1;k2;k3;k4];
trial_idx_m=[];
for i=1:4
    trial_idx_m=[trial_idx_m;ones(length(localmax_all_dat{1,i})/2,1);ones(length(localmax_all_dat{1,i})/2,1)*2];
end
type_idx_m=[ones(length(k1),1)*3;ones(length(k2),1)*4;ones(length(k3),1)*1;ones(length(k4),1)*2];
m_idx_m=[m1;m2;m3;m4];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show field example
field_all_dat={};
field_all_dat{1,1}=merge_cell(field_example(:,AD_idx==1)'); % Ntg old
field_all_dat{1,2}=merge_cell(field_example(:,AD_idx==2)');
field_all_dat{1,3}=merge_cell(field_example(:,AD_idx==3)');
field_all_dat{1,4}=merge_cell(field_example(:,AD_idx==4)');

k1=acorr_all_dat{1,1};
k2=acorr_all_dat{1,2};
idx=find(k2>1000)
k2(idx)=[]; % artifact, only 1 val >1000 while others <300
k3=acorr_all_dat{1,3};
k4=acorr_all_dat{1,4};

for i=1:length(k2)
    if(k3(i)>k4(i)&&k4(i)>k2(i)&&k2(i)>k1(i))
        idxx=[idxx,i];
    end
end

for k=1:length(idxx)
    figure;
    cell_list=[84,25,21]; % check area distribution by ploting the dat
    for i=1:4

        [bx1,by1]=boundary_from_binImg(field_all_dat{1,i}{idxx(k)}>=max(field_all_dat{1,i}{idxx(k)}(:))*0.5);
    %     [bx2,by2]=boundary_from_binImg(field_all_dat{1,i}{cell_list(2)}>=max(field_all_dat{1,i}{cell_list(2)}(:))*0.5);
    %     [bx3,by3]=boundary_from_binImg(field_all_dat{1,i}{cell_list(3)}>=max(field_all_dat{1,i}{cell_list(3)}(:))*0.5);

        subplot(3,4,i);imagesc(field_all_dat{1,i}{idxx(k)});hold on; plot(bx1,by1,'--','color','g'); % Ntg old
    %     subplot(3,4,i+4);imagesc(field_all_dat{1,i}{cell_list(2)});hold on; plot(bx2,by2,'--','color','g');% AD old
    %     subplot(3,4,i+8);imagesc(field_all_dat{1,i}{cell_list(3)});hold on; plot(bx3,by3,'--','color','g');% Ntg Y
    end
    colormap(jet)
    saveas(gcf,['D:\Xiaoxiao_3xtg_data\Rounds 14\figS3 panels\OF_autocorr_example\',num2str(k),'.tif']);
    saveas(gcf,['D:\Xiaoxiao_3xtg_data\Rounds 14\figS3 panels\OF_autocorr_example\',num2str(k),'.eps'],'epsc');
    close
end

% ATTN: use[4, 14, 43]
%% LT activate bin fractions
load('F:\AD_3xtg_linear_track_013121\dat_info.mat');
load('F:\AD_3xtg_linear_track_013121\all_neuron_behav_laps_042221.mat');

AD_idx=[4 2 4 2 4 2 4 3 3 1 3 3 4 4 3 1 3 1 4 4 2 2 2 1 1];
folderName=foldername;

bad_mice=[];
tk_range=[1:length(folderName)];
tk_range(bad_mice)=[];

binsize=15;

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

for tk=1:25
    trials_selection(tk,:)=[1:4];
end

for i=1:size(all_neuron_simp_laps1,1)
    for j=1:size(all_neuron_simp_laps1,2)
        try
            all_neuron_num(i,j)=size(all_neuron_simp_laps1{i,j}{1,1}.C,1);
        catch
            all_neuron_num(i,j)=size(all_neuron_simp_laps1{i,j}{1,2}.C,1);
        end
    end
end

fraction_all={};

for tk=tk_range
    fraction_all{tk}=[];
    for j=trials_selection(tk,:)
        val=[];
        [fr_all_dir,ct_all_dir,ctime_all_dir] = calculatinglinearTrackRateMaps_laps(all_neuron_simp_laps1{tk,j},all_behav_laps1{tk,j},binsize,'S',all_neuron_num(tk,j));
        
        fr_all_1=fr_all_dir{1}{1}; % 
        ctime_all_1=ctime_all_dir{1}{1}; %
        
        fr_all_2=fr_all_dir{2}{1}; % 
        ctime_all_2=ctime_all_dir{2}{1}; %
        
        val1=[];
        val2=[];
        for k=1:length(fr_all_1)   
            val1(k,1) = sum(sum(fr_all_1{k},1)>0)/sum(sum(ctime_all_1,1)>0);
            val2(k,1) = sum(sum(fr_all_2{k},1)>0)/sum(sum(ctime_all_2,1)>0);
        end
        fraction_all{tk}=[fraction_all{tk},[val1;val2]];
    end
    fraction_all{tk}=[mean(fraction_all{tk}(:,1:3),2);fraction_all{tk}(:,4)];
end

fraction_all_dat={};
fraction_all_dat{1,1}=cell2mat(fraction_all(:,AD_idx==3)');
fraction_all_dat{1,2}=cell2mat(fraction_all(:,AD_idx==4)');
fraction_all_dat{1,3}=cell2mat(fraction_all(:,AD_idx==1)');
fraction_all_dat{1,4}=cell2mat(fraction_all(:,AD_idx==2)');

fraction_all_dat1=fill_nan_to_cellmat(fraction_all_dat); %T first four column: run

%% LT auto correlation stats
autocorr_LT_area_all={};
autocorr_LT_rad_all={};
autocorr_LT_localmax_all={};
autocorr_LT_field_exp={};
midx={};
for tk=tk_range
    autocorr_LT_area_all{tk}=[];
    autocorr_LT_rad_all{tk}=[];
    autocorr_LT_localmax_all{tk}=[];
    autocorr_LT_field_exp{tk}={};
    field_temp={};
    for j=trials_selection(tk,:)
        val=[];
        [fr_all_dir,ct_all_dir,ctime_all_dir] = calculatinglinearTrackRateMaps_laps(all_neuron_simp_laps1{tk,j},all_behav_laps1{tk,j},binsize,'S',all_neuron_num(tk,j));
        
        fr_all_1=fr_all_dir{1}{1}; % 
        ctime_all_1=ctime_all_dir{1}{1}; %
        
        fr_all_2=fr_all_dir{2}{1}; % 
        ctime_all_2=ctime_all_dir{2}{1}; %
        
        val1=[];
        val2=[];
        rad1=[];
        rad2=[];
        localmax1=[];
        localmax2=[];
        
        midx_temp=[];
        for k=1:length(fr_all_1)   
%             if size(fr_all_1{k},1)>size(fr_all_1{k},2) % vertical
%                 fr_all_1{k}=nanmean(fr_all_1{k},2);
%                 fr_all_2{k}=nanmean(fr_all_2{k},2);
%             else
%                 fr_all_1{k}=nanmean(fr_all_1{k},1);
%                 fr_all_2{k}=nanmean(fr_all_2{k},1);
%             end
            [rad1(k,1),~,~,val1(k,1),localmax1(k,1),field_temp{k,j}]=firing_field_radius_autocorr(fr_all_1{k},1);
            [rad1(k,2),~,~,val2(k,1),localmax2(k,1)]=firing_field_radius_autocorr(fr_all_2{k},1);
%             [radius(k,1)]=firing_field_radius_autocorr(fr_all{tk,j}{k},1);  
            midx_temp(k,1)=tk;
        end
        autocorr_LT_area_all{tk}=[autocorr_LT_area_all{tk},[val1;val2]];
        autocorr_LT_rad_all{tk}=[autocorr_LT_rad_all{tk},[rad1;rad2]];
        autocorr_LT_localmax_all{tk}=[autocorr_LT_localmax_all{tk},[localmax1;localmax2]];
    end
    autocorr_LT_area_all{tk}=[nanmean(autocorr_LT_area_all{tk}(:,1:3),2);autocorr_LT_area_all{tk}(:,4)];
    autocorr_LT_rad_all{tk}=[nanmean(autocorr_LT_rad_all{tk}(:,1:3),2);autocorr_LT_rad_all{tk}(:,4)];
    autocorr_LT_localmax_all{tk}=[nanmean(autocorr_LT_localmax_all{tk}(:,1:3),2);autocorr_LT_localmax_all{tk}(:,4)];
    autocorr_LT_field_exp{tk}=field_temp(:,2);
    midx{tk}=[midx_temp;midx_temp;midx_temp;midx_temp];
end

autocorr_LT_area_all_dat={};
autocorr_LT_area_all_dat{1,1}=cell2mat(autocorr_LT_area_all(:,AD_idx==1)');
autocorr_LT_area_all_dat{1,2}=cell2mat(autocorr_LT_area_all(:,AD_idx==2)');
autocorr_LT_area_all_dat{1,3}=cell2mat(autocorr_LT_area_all(:,AD_idx==3)');
autocorr_LT_area_all_dat{1,4}=cell2mat(autocorr_LT_area_all(:,AD_idx==4)');

autocorr_LT_area_all_dat1=fill_nan_to_cellmat(autocorr_LT_area_all_dat); %T first four column: run

k1=autocorr_LT_area_all_dat{1,1};
k2=autocorr_LT_area_all_dat{1,2};
k3=autocorr_LT_area_all_dat{1,3};
k4=autocorr_LT_area_all_dat{1,4};

k1(k1==37)=[];
k2(k2==37)=[];
k3(k3==37)=[];
k4(k4==37)=[];
k1(k1==17)=[];
k2(k2==17)=[];
k3(k3==17)=[];
k4(k4==17)=[];

colorClusters_all=distinguishable_colors(10);

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[ 0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);
p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',30)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',30)

ranksum(k1,k2)
ranksum(k3,k4)
%% LT auto correlation rad
autocorr_LT_rad_all_dat{1,1}=cell2mat(autocorr_LT_rad_all(:,AD_idx==3)');
autocorr_LT_rad_all_dat{1,2}=cell2mat(autocorr_LT_rad_all(:,AD_idx==4)');
autocorr_LT_rad_all_dat{1,3}=cell2mat(autocorr_LT_rad_all(:,AD_idx==1)');
autocorr_LT_rad_all_dat{1,4}=cell2mat(autocorr_LT_rad_all(:,AD_idx==2)');

k1=cell2mat(autocorr_LT_rad_all(:,AD_idx==1)');
k2=cell2mat(autocorr_LT_rad_all(:,AD_idx==2)');
k3=cell2mat(autocorr_LT_rad_all(:,AD_idx==3)');
k4=cell2mat(autocorr_LT_rad_all(:,AD_idx==4)');

k1=floor(k1*1000)/1000;
k2=floor(k2*1000)/1000;
k3=floor(k3*1000)/1000;
k4=floor(k4*1000)/1000;
% 
% 
% k1(k1==[2.326])=[];
% k2(k2==[2.326])=[];
% k3(k3==[2.326])=[];
% k4(k4==[2.326])=[];
% 
% k1(k1==[2.705])=[];
% k2(k2==[2.705])=[];
% k3(k3==[2.705])=[];
% k4(k4==[2.705])=[];
% 
% 
% k1(k1==[3.6996])=[];
% k2(k2==[3.6996])=[];
% k1(k1==[3.8679])=[];
% k2(k2==[3.8679])=[];
% k3(k3==[4.0291])=[];
% k4(k4==[4.0291])=[];
% 
% k3(k3==4.029)=[];

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);
set(h3,'color',[ 0 1 0]);
set(h4,'color',[222/255,139/255,249/255]);
p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',30)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',30)

ranksum(k1,k2)
ranksum(k3,k4)
%% show LT field example
field_all_dat={};
field_all_dat{1,1}=merge_cell(autocorr_LT_field_exp(:,AD_idx==1)'); % Ntg old
field_all_dat{1,2}=merge_cell(autocorr_LT_field_exp(:,AD_idx==2)');
field_all_dat{1,3}=merge_cell(autocorr_LT_field_exp(:,AD_idx==3)');
field_all_dat{1,4}=merge_cell(autocorr_LT_field_exp(:,AD_idx==4)');

k1=autocorr_LT_area_all_dat{1,1};
k2=autocorr_LT_area_all_dat{1,2};
k3=autocorr_LT_area_all_dat{1,3};
k4=autocorr_LT_area_all_dat{1,4};

k1(k1==37)=[];
k2(k2==37)=[];
k3(k3==37)=[];
k4(k4==37)=[];
k1(k1==17)=[];
k2(k2==17)=[];
k3(k3==17)=[];
k4(k4==17)=[];

idxx=[];
for i=1:length(k2)
    if(k3(i)>k4(i)&&k1(i)>k3(i)&&k2(i)>k1(i))
        idxx=[idxx,i];
    end
end

for k=1:length(idxx)
    figure;
    try
    for i=1:4

        [bx1,by1]=boundary_from_binImg(field_all_dat{1,i}{idxx(k)}>=max(field_all_dat{1,i}{idxx(k)}(:))*0.5);

        subplot(3,4,i);imagesc(field_all_dat{1,i}{idxx(k)});hold on; plot(bx1,by1,'--','color','g'); % Ntg old
    end
    colormap(jet)
    saveas(gcf,['D:\Xiaoxiao_3xtg_data\Rounds 14\figS3 panels\LT_autocorr_example\',num2str(k),'.tif']);
    saveas(gcf,['D:\Xiaoxiao_3xtg_data\Rounds 14\figS3 panels\LT_autocorr_example\',num2str(k),'.eps'],'epsc');
    catch
    end
    close
end
