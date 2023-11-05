%% suppl 4: firing rate comparison

% Ntg young
load('D:\Xiaoxiao_3xtg_data\Rounds 14\figS4 panels\OF_FR.mat');
OF_fr_dat=cell2mat(fr_amp_dat1(1:2,:));
load('D:\Xiaoxiao_3xtg_data\Rounds 14\figS4 panels\LT_FR.mat');
LT_fr_dat=cell2mat(fr_amp_dat1(1:2,:));

colorClusters_all=distinguishable_colors(10);
% data order column: 3,4,1,2 (Ny, N

% Ntg young
k1=OF_fr_dat(:,1);
k2=LT_fr_dat(:,1);

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)

% ad young
k1=OF_fr_dat(:,2);
k2=LT_fr_dat(:,2);

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)

% Ntg old
k1=OF_fr_dat(:,3);
k2=LT_fr_dat(:,3);

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)

% AD old
k1=OF_fr_dat(:,4);
k2=LT_fr_dat(:,4);

h1=cdfplot(k1);hold on
h2=cdfplot(k2);
set(h1,'color',[0 0 1]);
set(h2,'color',[1 0 0]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',30)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',30)
