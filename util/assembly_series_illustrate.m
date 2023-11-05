function assembly_series_illustrate(timeSeries,group)
figure;
color_clust=distinguishable_colors(10);

if size(timeSeries,1)<=20 % 20 cells are too less, probably this is assembly series already
    for i=1:size(timeSeries,1)
        subplot(size(timeSeries,1),1,i);
        plot(timeSeries(i,:));
    end
else
    ugp=unique(group);
    ugp(ugp==-1)=[];
    
    for i=1:length(ugp)
        subplot(length(ugp),1,i);
        plot(mean(timeSeries(group==ugp(i),:),1),'color',color_clust(i,:));
    end
end

