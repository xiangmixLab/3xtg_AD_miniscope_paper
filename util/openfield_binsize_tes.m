p_summary={};

binsize_list=[5:5:50]
for binsize=1:length(binsize_list)
    try
    all_infoscore={};
    tic;
    for i=1:32
        for j=1:8
            [ii] = infoscore_check_021021(all_neuron_simp{i,j},all_behav{i,j}.position,all_behav{i,j}.time,'S',0,binsize_list(binsize),5);  
            all_infoscore{i,j}=ii(:,2);
            toc;
        end
    end

    ai=avg_infoscore(all_infoscore,trials_selection);

    info_per_type={ai(AD_idx==1),ai(AD_idx==2),ai(AD_idx==3),ai(AD_idx==4)};

    p_summary{binsize,1}=[ranksum(info_per_type{1},info_per_type{2}),ranksum(info_per_type{3},info_per_type{4}),ranksum(info_per_type{3},info_per_type{1}),ranksum(info_per_type{2},info_per_type{4})]

    catch
        disp([num2str(binsize),' not work'])
    end
end
