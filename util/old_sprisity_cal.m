spr_circle={};
spr_square={};

for tk=tk_range
    
    t1=[];
    t2=[];
    
    for i=trials_selection(tk,:)
        t=all_sparsity{tk,i};
        
        if ~isempty(findstr(cond_labels{i},'Circle'))&&~isempty(all_sparsity{tk,i})
                       
            t1=[t1;t];
        end
        if ~isempty(findstr(cond_labels{i},'Square'))&&~isempty(all_sparsity{tk,i})
                        
            t2=[t2;t];

        end        
    end
    spr_circle{tk,AD_idx(tk)}=t1;
    spr_square{tk,AD_idx(tk)}=t2;

end

colorClusters_all=distinguishable_colors(10);
% figure

% part 1: mice type discrimination
% subplot(121)
figure;
k1=[cell2mat(spr_circle(:,1));cell2mat(spr_square(:,1))];
k2=[cell2mat(spr_circle(:,2));cell2mat(spr_square(:,2))]+0.01;
k3=[cell2mat(spr_circle(:,3));cell2mat(spr_square(:,3))];
k4=[cell2mat(spr_circle(:,4));cell2mat(spr_square(:,4))];

all_size_inuse={k1,k2,k3,k4};
figure;
h1=cdfplot(k1);
hold on;
h2=cdfplot(k2);
h3=cdfplot(k3);
h4=cdfplot(k4);
% xlim([0 0.8])
set(h1,'color',colorClusters_all(1,:));
set(h2,'color',colorClusters_all(2,:));
set(h3,'color',colorClusters_all(3,:));
set(h4,'color',[222/255,139/255,249/255]);

p1=infer_cdf_loc(k1,nanmedian(k1));
p2=infer_cdf_loc(k2,nanmedian(k2));
p3=infer_cdf_loc(k3,nanmedian(k3));
p4=infer_cdf_loc(k4,nanmedian(k4));

plot(p1(1),p1(2),'.','color',colorClusters_all(1,:),'MarkerSize',20)
plot(p2(1),p2(2),'.','color',colorClusters_all(2,:),'MarkerSize',20)
plot(p3(1),p3(2),'.','color',colorClusters_all(3,:),'MarkerSize',20)
plot(p4(1),p4(2),'.','color',colorClusters_all(5,:),'MarkerSize',20)

[~,p]=kstest2(k1,k2)
[~,p]=kstest2(k3,k4)
[~,p]=kstest2(k1,k3)
[~,p]=kstest2(k2,k4)

ranksum(k1,k2)
ranksum(k3,k4)
ranksum(k1,k3)
ranksum(k2,k4)

nanmedian(k3)
nanmedian(k4)
nanmedian(k1)
nanmedian(k2)