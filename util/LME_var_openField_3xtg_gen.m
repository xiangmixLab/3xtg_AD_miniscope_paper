function [var,trial_idx_m,type_idx_m,m_idx_m]=LME_var_openField_3xtg_gen(dat,per_mice_dat,pc_idx,trial_selection,AD_idx)

type_idx={};
trial_idx={};
m_idx={};

infotype=2;


for i=1:size(per_mice_dat,1)
    
    if size(per_mice_dat,2)>1
        for j=1:size(trial_selection,2)
            type_idx{i,j}=ones(size(per_mice_dat{i,trial_selection(i,j)},1),1)*AD_idx(i);
            trial_idx{i,j}=ones(size(per_mice_dat{i,trial_selection(i,j)},1),1)*j;
            m_idx{i,j}=ones(size(per_mice_dat{i,trial_selection(i,j)},1),1)*i;
        end
    else
        for j=1:1
            type_idx{i,j}=ones(size(per_mice_dat{i,1},1),1)*AD_idx(i);
            trial_idx{i,j}=ones(size(per_mice_dat{i,1},1),1)*j;
            m_idx{i,j}=ones(size(per_mice_dat{i,1},1),1)*i;
        end
    end

end

if ~isempty(pc_idx)
    type_idx=infoscore_pc_trim(type_idx,pc_idx,infotype);
    trial_idx=infoscore_pc_trim(trial_idx,pc_idx,infotype);
    m_idx=infoscore_pc_trim(m_idx,pc_idx,infotype);
end

type_idx_all={cell2mat(reshape(type_idx(AD_idx==1,:),[],1)),cell2mat(reshape(type_idx(AD_idx==2,:),[],1)),cell2mat(reshape(type_idx(AD_idx==3,:),[],1)),cell2mat(reshape(type_idx(AD_idx==4,:),[],1))};

trial_idx_all={cell2mat(reshape(trial_idx(AD_idx==1,:),[],1)),cell2mat(reshape(trial_idx(AD_idx==2,:),[],1)),cell2mat(reshape(trial_idx(AD_idx==3,:),[],1)),cell2mat(reshape(trial_idx(AD_idx==4,:),[],1))};

m_idx_all={cell2mat(reshape(m_idx(AD_idx==1,:),[],1)),cell2mat(reshape(m_idx(AD_idx==2,:),[],1)),cell2mat(reshape(m_idx(AD_idx==3,:),[],1)),cell2mat(reshape(m_idx(AD_idx==4,:),[],1))};

k1=dat{1};
k2=dat{2};
k3=dat{3};
k4=dat{4};

var=[k1;k2;k3;k4];
trial_idx_m=[];
type_idx_m=[];
m_idx_m=[];

for i=1:4
 
    trial_idx_m=[trial_idx_m;trial_idx_all{i}(:,1)];
    type_idx_m=[type_idx_m;type_idx_all{i}(:,1)];
    m_idx_m=[m_idx_m;m_idx_all{i}(:,1)];

end

nanvar_idx=isnan(var);
var(nanvar_idx)=[];
trial_idx_m(nanvar_idx)=[];
type_idx_m(nanvar_idx)=[];
m_idx_m(nanvar_idx)=[];
