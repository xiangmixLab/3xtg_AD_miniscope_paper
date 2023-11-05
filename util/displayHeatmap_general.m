function ax=displayHeatmap_general(CM1,CM2,perm,colorClusters,optimalK,group,flipSign)

ax=figure;
set(gcf,'outerposition',get(0,'screensize'));

subplot(1,10,[5 10])

CM1(logical(eye(size(CM1)))) = 1;
if flipSign==0
    imagesc(flip(CM1(perm,perm)));
else
    imagesc(CM1(perm,perm));
end
% axis square
set(gca,'Xtick',[]);set(gca,'Ytick',[])
colormap(jet)
caxis([min(CM1(:)) max(CM1(:))])
c = colorbar;
c.Location = 'EastOutside';
c.Label.String = 'Simimarity';
c.Label.FontSize = 8;%c.Label.FontWeight = 'bold'
c.FontSize = 8;
title('CorrelationHeatMap');


subplot(1,10,[1 2])
% optimalK = length(colorClusters);
% Z = linkage(CM2,'complete');
% group = cluster(Z,'maxclust',optimalK);
% 
% H = dendrogram(Z,0,'Orientation','left','Labels',[],'reorder',perm);  %[H,T,PERM] = dendrogram(Z,0,'colorthreshold',0.35); %dendrogram plot with colorthreshold set to 0.35
% set(gca,'fontsize',8);
% set(H,'color','k')
% ylim([0.5 length(group)+0.5])
% set(gca,'Xtick',[]);set(gca,'Ytick',[]);


% 
subplot(1,10,[3 4])
groupUni = unique(group(perm),'stable');
ybar = [];
for i = 1:1:optimalK
    ybar = [ybar,nnz(group == groupUni(i))];
end
H = bar([ybar;nan(1,length(ybar))],'stacked','BarWidth',0.0004,'ShowBaseLine','off','LineWidth',0.01);
k = 0;
for i = 1:1:optimalK
    k = k+1;
    set(H(i),'facecolor',colorClusters(groupUni(k),:))
end
axis tight
xlim([1 1.001])
set(gca,'xtick',1)
axis off
if flipSign==0
else
    axis ij;
end

colormap(jet)
caxis([0,0.8*max(CM1(:))])