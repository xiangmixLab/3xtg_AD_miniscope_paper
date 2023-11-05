function [traj_leng,behavpos_dis_list]=trajectory_length(behavpos)

behavpos_dis_list=zeros(size(behavpos,1)-1,1);
for i=1:size(behavpos,1)-1
    behavpos_dis_list(i)=sum((behavpos(i+1,:)-behavpos(i,:)).^2).^0.5;
end

traj_leng=nansum(behavpos_dis_list);