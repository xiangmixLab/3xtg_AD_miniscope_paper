
function [all_dis_corr_circle_all_intra,all_dis_corr_circle_all_inter,all_dis_corr_circle_all_all,midx,f,gof]=pairwise_dis_tempCorr_031722(foldername,group1,session)

all_dis_circle_intra={};
all_corr_circle_intra={};
all_dis_circle_inter={};
all_corr_circle_inter={};
all_dis_circle_all={};
all_corr_circle_all={};


for k=1:length(foldername)
    cd(foldername{k});
    load('further_processed_neuron_extraction_final_result.mat');
    load('neuronIndividuals_new.mat');
%     load([cluster_file_name{k}])
    group=group1{k,session};              
    uni_group=unique(group);
    uni_group(uni_group==-1)=[];
    
    neuronIndividuals_new{session}.centroid=neuron_centroid_calculation(neuronIndividuals_new{session},neuronIndividuals_new{session}.imageSize);
    dis_all=squareform(pdist(neuronIndividuals_new{session}.centroid))*2; % *2 is added 042320, in light of the resampling of the frame make all centroid distance only half of what it should be
    corr_all=squareform(1-pdist(neuronIndividuals_new{session}.C,'correlation'));
    
    % all
    all_corr_circle_all{k}=corr_all;
    all_dis_circle_all{k}=dis_all;
    % in group
    tic;
    for k1=1:length(uni_group)
        %distance
        t=triu(dis_all(group==k1,group==k1),1);
        mask = tril(true(size(corr_all(group==k1,group==k1))),1);
        all_dis_circle_intra{k,k1}=t(~mask);
        t=triu(dis_all(group==k1,group~=k1),1);
        mask = tril(true(size(corr_all(group==k1,group~=k1))),1);
        all_dis_circle_inter{k,k1}=t(~mask);
        %correlation
        t=triu(corr_all(group==k1,group==k1),1);
        mask = tril(true(size(corr_all(group==k1,group==k1))),1);
        all_corr_circle_intra{k,k1}=t(~mask);
        t=triu(corr_all(group==k1,group~=k1),1);
        mask = tril(true(size(corr_all(group==k1,group~=k1))),1);
        all_corr_circle_inter{k,k1}=t(~mask);   
        toc;
    end
end

all_dis_corr_circle_all_intra=[];
all_dis_corr_circle_all_inter=[];
all_dis_corr_circle_all_all=[];
ctt1=1;
ctt2=1;
ctt3=1;

midx={};
midx_intra=[];
midx_inter=[];
midx_all=[];

dis_thresh=20;

for k=1:length(foldername)
    cd(foldername{k});
    group=group1{k,session};              
    uni_group=unique(group);
    uni_group(uni_group==-1)=[];  
    
    for i=1:length(uni_group)
        
        dis_non_close_intra=all_dis_circle_intra{k,i}>dis_thresh; % 35um distance
        all_dis_corr_circle_all_intra(ctt1:ctt1+length(all_dis_circle_intra{k,i}(dis_non_close_intra))-1,1)=all_dis_circle_intra{k,i}(dis_non_close_intra);
        all_dis_corr_circle_all_intra(ctt1:ctt1+length(all_dis_circle_intra{k,i}(dis_non_close_intra))-1,2)=all_corr_circle_intra{k,i}(dis_non_close_intra);
        midx_intra(ctt1:ctt1+length(all_dis_circle_intra{k,i}(dis_non_close_intra))-1)=k;
        ctt1=ctt1+length(all_dis_circle_intra{k,i}(dis_non_close_intra));
        
        dis_non_close_inter=all_dis_circle_inter{k,i}>dis_thresh; % 35um distance
        all_dis_corr_circle_all_inter(ctt2:ctt2+length(all_dis_circle_inter{k,i}(dis_non_close_inter))-1,1)=all_dis_circle_inter{k,i}(dis_non_close_inter);
        all_dis_corr_circle_all_inter(ctt2:ctt2+length(all_dis_circle_inter{k,i}(dis_non_close_inter))-1,2)=all_corr_circle_inter{k,i}(dis_non_close_inter);
        midx_inter(ctt2:ctt2+length(all_dis_circle_inter{k,i}(dis_non_close_inter))-1)=k;
        ctt2=ctt2+length(all_dis_circle_inter{k,i}(dis_non_close_inter));
    end
    
    dis_non_close_all=all_dis_circle_all{k}>dis_thresh; % 35um distance
    all_dis_corr_circle_all_all(ctt3:ctt3+length(all_dis_circle_all{k}(dis_non_close_all))-1,1)=all_dis_circle_all{k}(dis_non_close_all);
    all_dis_corr_circle_all_all(ctt3:ctt3+length(all_dis_circle_all{k}(dis_non_close_all))-1,2)=all_corr_circle_all{k}(dis_non_close_all);
    midx_all(ctt3:ctt3+length(all_dis_circle_all{k}(dis_non_close_all))-1)=k;
    ctt3=ctt3+length(all_dis_circle_all{k}(dis_non_close_all));

end

midx={midx_intra,midx_inter,midx_all};

non_nan_idx_intra=logical(~isnan(all_dis_corr_circle_all_intra(1:end,1)).*~isnan(all_dis_corr_circle_all_intra(1:end,2)));
non_nan_idx_inter=logical(~isnan(all_dis_corr_circle_all_inter(1:end,1)).*~isnan(all_dis_corr_circle_all_inter(1:end,2)));
non_nan_idx_all=logical(~isnan(all_dis_corr_circle_all_all(1:end,1)).*~isnan(all_dis_corr_circle_all_all(1:end,2)));

[f1,gof1]=fit(all_dis_corr_circle_all_intra(non_nan_idx_intra,1),all_dis_corr_circle_all_intra(non_nan_idx_intra,2),'power2'); % 031722 flipped x and y 
[f2,gof2]=fit(all_dis_corr_circle_all_inter(non_nan_idx_inter,1),all_dis_corr_circle_all_inter(non_nan_idx_inter,2),'power2');
[f3,gof3]=fit(all_dis_corr_circle_all_all(non_nan_idx_all,1),all_dis_corr_circle_all_all(non_nan_idx_all,2),'power2');

figure;
% plot(all_dis_corr_circle_all_intra(1:10:end,1),all_dis_corr_circle_all_intra(1:10:end,2),'.','MarkerSize',15,'color','b');
hold on;
% plot(all_dis_corr_circle_all_inter(1:50:end,1),all_dis_corr_circle_all_inter(1:50:end,2),'.','MarkerSize',15,'color','r');
h2=plot(f2,all_dis_corr_circle_all_inter(1:200:end,1),all_dis_corr_circle_all_inter(1:200:end,2));
h1=plot(f1,all_dis_corr_circle_all_intra(1:50:end,1),all_dis_corr_circle_all_intra(1:50:end,2));
h3=plot(f3);

set([h1],'color',[126,47,142]/255,'MarkerSize',10)
set([h2],'color',[119,172,48]/255,'MarkerSize',10)
set([h3],'color','c','MarkerSize',10)
set(gca,'ylim',[-0.2 1])
% [b1,~,~,~,stats1]=regress(all_dis_corr_linear_all(1:end,2),[all_dis_corr_linear_all(1:end,1),ones(length(all_dis_corr_linear_all(1:end,1)),1)]);
% plot(all_dis_corr_linear_all(1:10000:end,1),b1(2)+b1(1)*all_dis_corr_linear_all(1:10000:end,1),'r');

% all_dis_corr_circle_all_intra(isnan(all_dis_corr_circle_all_intra(:,1)),:)=[];
% all_dis_corr_circle_all_intra(isnan(all_dis_corr_circle_all_intra(:,2)),:)=[];
% all_dis_corr_circle_all_inter(isnan(all_dis_corr_circle_all_inter(:,1)),:)=[];
% all_dis_corr_circle_all_inter(isnan(all_dis_corr_circle_all_inter(:,2)),:)=[];
f={f1,f2,f3};

dat=all_dis_corr_circle_all_all(non_nan_idx_all,:);
dat=dat(1:100:end,:);
[rho,pval] = corr(dat(:,1),dat(:,2),'Type','Spearman');
gof={gof1,gof2,gof3,[rho,pval]};