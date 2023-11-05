function [var,trial_idx_m,type_idx_m,m_idx_m,dir_idx_m]=LME_var_linearTrack_3xtg_gen(dat,per_mice_dat_dir1,per_mice_dat_dir2,pc_idx_dir1,pc_idx_dir2,trialNum,AD_idx)

infotype=2;
type_dir1_idx={};
trial_dir1_idx={};
m_dir1_idx={};

type_dir2_idx={};
trial_dir2_idx={};
m_dir2_idx={};

dir1_idx={};
dir2_idx={};

for i=1:size(per_mice_dat_dir1,1)
    type_dir1_idx{i,1}=[];
    trial_dir1_idx{i,1}=[];
    m_dir1_idx{i,1}=[];
    dir1_idx{i,1}=[];

    type_dir2_idx{i,1}=[];
    trial_dir2_idx{i,1}=[];
    m_dir2_idx{i,1}=[];
    dir2_idx{i,1}=[];

    for j=1:trialNum
        type_dir1_tmp=ones(size(per_mice_dat_dir1{i},1),1)*AD_idx(i);
        trial_dir1_tmp=ones(size(per_mice_dat_dir1{i},1),1)*j;
        m_dir1_tmp=ones(size(per_mice_dat_dir1{i},1),1)*i;
        dir1_tmp=ones(size(per_mice_dat_dir1{i},1),1);
        
        type_dir2_tmp=ones(size(per_mice_dat_dir2{i},1),1)*AD_idx(i);
        trial_dir2_tmp=ones(size(per_mice_dat_dir2{i},1),1)*j;
        m_dir2_tmp=ones(size(per_mice_dat_dir2{i},1),1)*i;
        dir2_tmp=ones(size(per_mice_dat_dir2{i},1),1);

        if ~isempty(pc_idx_dir1) && ~isempty(pc_idx_dir2)
            if iscell(pc_idx_dir1{i,j})
                type_dir1_tmp=type_dir1_tmp(pc_idx_dir1{i,j}{infotype});
                trial_dir1_tmp=trial_dir1_tmp(pc_idx_dir1{i,j}{infotype});
                m_dir1_tmp=m_dir1_tmp(pc_idx_dir1{i,j}{infotype});
                dir1_tmp=dir1_tmp(pc_idx_dir1{i,j}{infotype});
    
                type_dir2_tmp=type_dir2_tmp(pc_idx_dir2{i,j}{infotype});
                trial_dir2_tmp=trial_dir2_tmp(pc_idx_dir2{i,j}{infotype});
                m_dir2_tmp=m_dir2_tmp(pc_idx_dir2{i,j}{infotype});
                dir2_tmp=dir2_tmp(pc_idx_dir2{i,j}{infotype});
            else
                type_dir1_tmp=type_dir1_tmp(pc_idx_dir1{i,j});
                trial_dir1_tmp=trial_dir1_tmp(pc_idx_dir1{i,j});
                m_dir1_tmp=m_dir1_tmp(pc_idx_dir1{i,j});
                dir1_tmp=dir1_tmp(pc_idx_dir1{i,j});
    
                type_dir2_tmp=type_dir2_tmp(pc_idx_dir2{i,j});
                trial_dir2_tmp=trial_dir2_tmp(pc_idx_dir2{i,j});
                m_dir2_tmp=m_dir2_tmp(pc_idx_dir2{i,j});
                dir2_tmp=dir2_tmp(pc_idx_dir2{i,j});
            end
        end

        type_dir1_idx{i,1}=[type_dir1_idx{i,1};type_dir1_tmp];
        trial_dir1_idx{i,1}=[trial_dir1_idx{i,1};trial_dir1_tmp];
        m_dir1_idx{i,1}=[m_dir1_idx{i,1};m_dir1_tmp];
        dir1_idx{i,1}=[dir1_idx{i,1};dir1_tmp];
    
        type_dir2_idx{i,1}=[type_dir2_idx{i,1};type_dir2_tmp];
        trial_dir2_idx{i,1}=[trial_dir2_idx{i,1};trial_dir2_tmp];
        m_dir2_idx{i,1}=[m_dir2_idx{i,1};m_dir2_tmp];
        dir2_idx{i,1}=[dir2_idx{i,1};dir2_tmp];
    end
end

type_dir1_idx_all={cell2mat(type_dir1_idx(AD_idx==1)),cell2mat(type_dir1_idx(AD_idx==2)),cell2mat(type_dir1_idx(AD_idx==3)),cell2mat(type_dir1_idx(AD_idx==4))};
type_dir2_idx_all={cell2mat(type_dir2_idx(AD_idx==1)),cell2mat(type_dir2_idx(AD_idx==2)),cell2mat(type_dir2_idx(AD_idx==3)),cell2mat(type_dir2_idx(AD_idx==4))};

trial_dir1_idx_all={cell2mat(trial_dir1_idx(AD_idx==1)),cell2mat(trial_dir1_idx(AD_idx==2)),cell2mat(trial_dir1_idx(AD_idx==3)),cell2mat(trial_dir1_idx(AD_idx==4))};
trial_dir2_idx_all={cell2mat(trial_dir2_idx(AD_idx==1)),cell2mat(trial_dir2_idx(AD_idx==2)),cell2mat(trial_dir2_idx(AD_idx==3)),cell2mat(trial_dir2_idx(AD_idx==4))};

m_dir1_idx_all={cell2mat(m_dir1_idx(AD_idx==1)),cell2mat(m_dir1_idx(AD_idx==2)),cell2mat(m_dir1_idx(AD_idx==3)),cell2mat(m_dir1_idx(AD_idx==4))};
m_dir2_idx_all={cell2mat(m_dir2_idx(AD_idx==1)),cell2mat(m_dir2_idx(AD_idx==2)),cell2mat(m_dir2_idx(AD_idx==3)),cell2mat(m_dir2_idx(AD_idx==4))};

dir1_idx_all={cell2mat(dir1_idx(AD_idx==1)),cell2mat(dir1_idx(AD_idx==2)),cell2mat(dir1_idx(AD_idx==3)),cell2mat(dir1_idx(AD_idx==4))};
dir2_idx_all={cell2mat(dir2_idx(AD_idx==1)),cell2mat(dir2_idx(AD_idx==2)),cell2mat(dir2_idx(AD_idx==3)),cell2mat(dir2_idx(AD_idx==4))};

k1=dat{1};
k2=dat{2};
k3=dat{3};
k4=dat{4};

var=[k1;k2;k3;k4];
trial_idx_m=[];
type_idx_m=[];
m_idx_m=[];
dir_idx_m=[];

for i=1:4

    trial_idx_m=[trial_idx_m;[trial_dir1_idx_all{i}(:,1);trial_dir2_idx_all{i}(:,1)]];
    type_idx_m=[type_idx_m;[type_dir1_idx_all{i}(:,1);type_dir2_idx_all{i}(:,1)]];
    m_idx_m=[m_idx_m;[m_dir1_idx_all{i}(:,1);m_dir2_idx_all{i}(:,1)]];
    dir_idx_m=[dir_idx_m;[dir1_idx_all{i}(:,1);dir2_idx_all{i}(:,1)]];

end

nanvar_idx=isnan(var);
var(nanvar_idx)=[];
trial_idx_m(nanvar_idx)=[];
type_idx_m(nanvar_idx)=[];
m_idx_m(nanvar_idx)=[];
dir_idx_m(nanvar_idx)=[];
