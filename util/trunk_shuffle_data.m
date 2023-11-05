function dat1=trunk_shuffle_data(dat,trunk_num)

interval=round(size(dat,2)/(trunk_num));
interval_list=[ones(1,floor(size(dat,2)/interval))*interval,size(dat,2)-floor(size(dat,2)/interval)*interval];
interval_list(interval_list==0)=[]; % the last one is 0 inperfact case

dat_trunk=mat2cell(dat',interval_list); %mat2cell(dat,[row_sizes_of_the_mat_in_current_splitted_cell],[col_sizes_of_the_mat_in_current_splitted_cell])

rand_order=randperm(size(dat_trunk,1));
dat_shuf=cell2mat(dat_trunk(rand_order));
dat1=dat_shuf';