function [intra_all,inter_all,avg_intra_all,avg_inter_all,intra_shuffle_all,avg_intra_shuffle_all,avg_intra_shuffle_all_by1000,intra_all_neuronBased,inter_all_neuronBased,intra_shuffle_all_neuronBased]=intra_inter_cluster_corr_dis(neuronIndividuals,group,session,sign)

shuffle=1000;

intra={};
inter={};
intra_shuffle={};

intra_neuron={};
inter_neuron={};
intra_shuffle_neuron={};

ugp=unique(group);
ugp(ugp==-1)=[];

switch sign
    case 'corr'
        dataC1=neuronIndividuals{session}.C;
        dataC1(isnan(dataC1(:)))=0;
        dataC1=zscore(dataC1,[],2);
        
        corr_all=squareform(1-pdist(dataC1,'correlation'));
    case 'dis'
        dataC1=neuron_centroid_calculation(neuronIndividuals{session},size(neuronIndividuals{session}.Cn));
        corr_all=squareform(pdist(dataC1))*2; % if using 240*376 option, the distance is cut in half      
    otherwise
        dataC1=neuronIndividuals{session}.C;
        dataC1(isnan(dataC1(:)))=0;
        dataC1=zscore(dataC1,[],2);
        
        corr_all=squareform(1-pdist(dataC1,'correlation'));
end

% cell index
cellIdx=corr_all*0;
for i=1:size(cellIdx,1)
    for j=i+1:size(cellIdx,2)
        cellIdx(i,j)=i;
    end
end

for i=1:length(ugp)
    t=triu(corr_all(group==i,group==i),1);
    cellIdxt=triu(cellIdx(group==i,group==i),1);
    mask = tril(true(size(corr_all(group==i,group==i))),1);
    intra{i}=t(~mask);
    intra_neuron{i}=cellIdxt(~mask);
    
    t=triu(corr_all(group==i,group~=i),1);
    cellIdxt=triu(cellIdx(group==i,group~=i),1);
    mask = tril(true(size(corr_all(group==i,group~=i))),1);
    inter{i}=t(~mask);
    inter_neuron{i}=cellIdxt(~mask);
end

for j=1:shuffle
    group_rand = group(randperm(length(group)));
    for i=1:length(ugp)
        t=triu(corr_all(group_rand==i,group_rand==i),1); 
        cellIdxt=triu(cellIdx(group_rand==i,group_rand==i),1);
        mask = tril(true(size(corr_all(group_rand==i,group_rand==i))),1);
        intra_shuffle{i,j}=t(~mask);
        intra_shuffle_neuron{i,j}=cellIdxt(~mask);
    end
end
    
intra_all=[];
ctt=1;
for i=1:length(ugp)
    intra_all(ctt:ctt+length(intra{i})-1)=intra{i};
    ctt=ctt+length(intra{i});
end

inter_all=[];
ctt=1;
for i=1:length(ugp)
    inter_all(ctt:ctt+length(inter{i})-1)=inter{i};
    ctt=ctt+length(inter{i});
end

intra_shuffle_all=[];
ctt=1;
for j=1:shuffle
    for i=1:length(ugp)
        intra_shuffle_all(ctt:ctt+length(intra_shuffle{i,j})-1)=intra_shuffle{i,j};
        ctt=ctt+length(intra_shuffle{i,j});
    end
end

avg_intra_shuffle_all_by1000=[];
for i=1:length(ugp)
    as_t=[];
    for j=1:shuffle
        as_t(:,j)=intra_shuffle{i,j};
    end
    avg_intra_shuffle_all_by1000=[avg_intra_shuffle_all_by1000;nanmean(as_t,2)];
end

intra_all(isnan(intra_all))=[];
inter_all(isnan(inter_all))=[];
intra_shuffle_all(isnan(intra_shuffle_all))=[];

intra_all=intra_all';
inter_all=inter_all';
intra_shuffle_all=intra_shuffle_all';

avg_intra_all=nanmean(intra_all);
avg_inter_all=nanmean(inter_all);
avg_intra_shuffle_all=nanmean(intra_shuffle_all);

% neuron based
intra_all_neuronBased=[];
ctt=1;
for i=1:length(ugp)
    intra_tmp=[];
    u_cell=unique(intra_neuron{i});
    for j=1:length(unique(u_cell))
        intra_tmp(j,1)=nanmean(intra{i}(intra_neuron{i}==u_cell(j)));
    end
    intra_all_neuronBased(ctt:ctt+length(intra_tmp)-1)=intra_tmp;
    ctt=ctt+length(intra_tmp);
end

inter_all_neuronBased=[];
ctt=1;
for i=1:length(ugp)
    inter_tmp=[];
    u_cell=unique(inter_neuron{i});
    for j=1:length(unique(u_cell))
        inter_tmp(j,1)=nanmean(inter{i}(inter_neuron{i}==u_cell(j)));
    end
    inter_all_neuronBased(ctt:ctt+length(inter_tmp)-1)=inter_tmp;
    ctt=ctt+length(inter_tmp);
end

intra_shuffle_all_neuronBased=[];
ctt=1;

for i=1:length(ugp)
    u_cell=unique(intra_shuffle_neuron{i,1});
    intra_shuf_tmp=zeros(length(unique(u_cell)),1);
    countt=zeros(length(unique(u_cell)),1);
    
    for j=1:shuffle

        for k=1:length(unique(u_cell))
            if ~isnan(nanmean(intra_shuffle{i,j}(intra_shuffle_neuron{i,j}==u_cell(k))))
                intra_shuf_tmp(k,1)=intra_shuf_tmp(k,1)+nanmean(intra_shuffle{i,j}(intra_shuffle_neuron{i,j}==u_cell(k)));
                countt(k,1)=countt(k,1)+1;
            end
        end
        
    end
    intra_shuf_tmp=intra_shuf_tmp./countt;
    intra_shuffle_all_neuronBased(ctt:ctt+length(intra_shuf_tmp)-1)=intra_shuf_tmp;
    ctt=ctt+length(intra_shuf_tmp);
end

intra_all_neuronBased(isnan(intra_all_neuronBased))=[];
inter_all_neuronBased(isnan(inter_all_neuronBased))=[];
intra_shuffle_all_neuronBased(isnan(intra_shuffle_all_neuronBased))=[];

intra_all_neuronBased=intra_all_neuronBased';
inter_all_neuronBased=inter_all_neuronBased';
intra_shuffle_all_neuronBased=intra_shuffle_all_neuronBased';

