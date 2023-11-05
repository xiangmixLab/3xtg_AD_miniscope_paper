function clustered_timeSeries2(ts,group2)

% dataC2 = [];
% for j = 1:length(unique(group2))
%     dataC2 = [dataC2;ts(group2 == j,:)];
% end
% 
% positionC2 = 0;
% for i =1:length(unique(group2))
%     positionC2(i+1) = positionC2(i)+sum(group2 == i);
% end

% normalize ts
ts=ts./max(ts,[],2)*20;
% imshow 
clust_inuse=unique(group2);
clust_inuse(clust_inuse==-1)=[];
figure;
colorClusters2=distinguishable_colors(10);
for j = 1:length(clust_inuse)

    subplot(length(clust_inuse),1,j);
    
    rangee=find(group2 == clust_inuse(j));
    rangee=rangee(1:2:end); % only show half
    ctt=1;
    for k=1:length(rangee)
        plot(ts(rangee(k),1:end)'+(ctt-1)*1,'-','color',colorClusters2(j,:));hold on;
        ctt=ctt+1;
    end
    axis off
end
 set(gcf,'renderer','painters');
