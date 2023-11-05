% this cluster overlap is still based on rand index, but follows modi et
% al. 2014.
% original rand index = (a+b)/nchoosek(n,2), where a is the number of pairs
% stay together in both and b is the # of pairs separate in both
% This overlap is just (a)/nchoosek(n,2). no separate pairs considered
% 
function [overlap_idx_avgClust,overlap_idx,overlap_idx_perClust]=new_cluster_overlap_latest(c1,c2)

C=Contingency(c1,c2);	%form contingency matrix

%% total pairs of elements matched elements between c1 and c2
together_pairs=[];
for i=1:size(C,1)
    for j=1:size(C,2)
        if C(i,j)>2
            together_pairs(i,j)=nchoosek(C(i,j),2);
        end
        if C(i,j)==2 && C(i,j)>0
            together_pairs(i,j)=1;
        end
        if C(i,j)==0
            together_pairs(i,j)=0; % if this item is 0, equals to Rand index's a
        end
    end
end

clust_pairs=[];
for i=1:size(C,1)
    if sum(C(i,:))>2
        clust_pairs(i,1)=nchoosek(sum(C(:,i)),2);
    end
    if sum(C(i,:))==2
        clust_pairs(i,1)=1;
    end
    if sum(C(i,:))<2
        clust_pairs(i,1)=1; % only one element in the cluster, set as 1 to avoid divide by 0
    end
end

overlap_idx_avgClust=mean(sum(together_pairs,1)./clust_pairs'); % overlap for each cluster, then take average
overlap_idx=sum(together_pairs(:))/nchoosek(length(c1),2); % overlap calculated by all overlap pairs / all possible pairs
overlap_idx_perClust=sum(together_pairs,1)./clust_pairs';


function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
    if Mem1(i)>0&&Mem2(i)>0
        Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
    end
end
