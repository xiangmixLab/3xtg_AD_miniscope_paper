function intra_inter_correlation_cdf(intra_all,inter_all,intra_shuffle_all)

colorClusters_all=colormap(lines);close

figure;

h1=cdfplot(cell2mat(intra_all'));
hold on;
h2=cdfplot(cell2mat(inter_all'));
h3=cdfplot(cell2mat(intra_shuffle_all'));
set(h1,'color',colorClusters_all(4,:));
set(h2,'color',colorClusters_all(5,:));
set(h3,'color',colorClusters_all(6,:));

p1=infer_cdf_loc(cell2mat(intra_all'),mean(cell2mat(intra_all')));
p2=infer_cdf_loc(cell2mat(inter_all'),mean(cell2mat(inter_all')));
p3=infer_cdf_loc(cell2mat(intra_shuffle_all'),mean(cell2mat(intra_shuffle_all')));

plot(p1(1),p1(2),'.','color',colorClusters_all(4,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(6,:),'MarkerSize',20)

disp(['intra-cluster correlation mean:', num2str(mean(cell2mat(intra_all')))])
disp(['intra-cluster correlation standard error:',num2str(sem(cell2mat(intra_all'),1))])
disp(['inter-cluster correlation avg:',num2str(mean(cell2mat(inter_all')))])
disp(['inter-cluster correlation avg:',num2str(sem(cell2mat(inter_all'),1))])
disp(['shuffled intra-cluster correlation avg:',num2str(mean(cell2mat(intra_shuffle_all')))])
disp(['shuffled intra-cluster correlation sem:',num2str(sem(cell2mat(intra_shuffle_all'),1))])

p1=kstest2_sample(cell2mat(intra_all'),cell2mat(inter_all'),0.01);
p2=kstest2_sample(cell2mat(intra_all'),cell2mat(intra_shuffle_all'),0.01);

disp(['intra-cluster vs inter-cluster:', num2str(p1)])
disp(['intra-cluster vs shuffled_intra_cluster:',num2str(p2)])

