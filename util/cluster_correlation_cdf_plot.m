function cluster_correlation_cdf_plot(data,colorss)
figure; hold on;
for i=1:length(data)
    
    if iscell(data{i})
        h1=cdfplot(cell2mat(data{i}'));
        p1=infer_cdf_loc(cell2mat(data{i}'),mean(cell2mat(data{i}')));
    else
        h1=cdfplot(data);
        p1=infer_cdf_loc(data{i},mean(data{i}));
    end
    set(h1,'color',colorss(i,:));
    plot(p1(1),p1(2),'.','color',colorss(i,:),'MarkerSize',20)
end
