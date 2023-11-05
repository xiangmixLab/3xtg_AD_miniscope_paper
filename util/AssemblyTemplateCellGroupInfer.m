function group=AssemblyTemplateCellGroupInfer(AssemblyTemplatesIndividual,time_projectionIndividual,nC,weightThreshPercentile)

% at_z=zscore(AssemblyTemplatesIndividual,[],1); % zscore to each row 
at_z=AssemblyTemplatesIndividual;
at_z_idx=[];
at_z_sort=[];
for i=1:size(at_z,2)
    [at_z_sort(:,i),at_z_idx(:,i)]=sort(at_z(:,i),'descend');
end

for i=1:size(at_z,2)
    at=at_z_sort(:,i);
    if ~isempty(weightThreshPercentile)
        at(at<quantile(at,weightThreshPercentile))=-1;
    end      
    at_z_sort(:,i)=at;
end
% use size(at_z,2) pointers and a map
% each pointer represents a group
% map record all visited indexes and their corresponding component,
% position

% if there's a duplicate index, compare the order of the previous record
% and current record, assign the index to the group in which the index is
% in more important position

% if there's a tie (duplicate in current round), compare the correlation
% between the component and the neuron 

map=[-1,-1,-1];

for i=1:size(at_z_sort,1)
    for j=1:size(at_z_sort,2)
        if at_z_sort(i,j)>0 % only examine above avg contribution
            if ~ismember(map(:,1),at_z_idx(i,j))
                map=[map;[at_z_idx(i,j),j,i]]; % index, curGroup, priority position
            else
                % index already been recorded
                previousRecord=map(map(:,1)==at_z_idx(i,j),2:end);

                if i==previousRecord(end)
                    % tie
                    curComponent=time_projectionIndividual(j,:);
                    prevComponent=time_projectionIndividual(previousRecord(1),:);
                    curTransient=nC(at_z_idx(i,j),:);
                    
%                     curTransient=resample(curTransient, length(curComponent),length(curTransient));
                    
                    corrCur=corrcoef(curTransient,curComponent);
                    corrPrev=corrcoef(curTransient,prevComponent);
                    corrCur=corrCur(2);
                    corrPrev=corrPrev(2);

                    if corrCur>corrPrev
                        map(map(:,1)==at_z_idx(i,j),:)=[at_z_idx(i,j),j,i]; 
                    end

                else
                    % i>previousRecord(end), keep previous record
                end

            end
        end
    end
end

group=map(2:end,1:2);


remaining_idx=setdiff([1:size(nC,1)],group(:,1));
group=[group;[remaining_idx',ones(length(remaining_idx),1)*-1]];
            
group=sortrows(group,1);
group=group(:,2);