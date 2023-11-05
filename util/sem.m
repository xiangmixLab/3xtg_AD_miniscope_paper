function sem_res=sem(dat,dim)

sem_res=nanstd(dat,[],dim)/(size(dat,dim))^0.5;