function preserve_idxs=exclude_neurons_by_temporal(dpath)

preserve_idxs={};

for i=1:length(dpath)
    load([dpath{i},'\','neuronIndividuals_new_tuncat.mat']);
    
    for j=1:1
        preserve_idxs{i,j}=zeros(size(neuronIndividuals_new_tuncat{j}.C,1),1);
        for k=1:size(neuronIndividuals_new_tuncat{j}.C,1)
            plot(neuronIndividuals_new_tuncat{j}.C(k,:));
            strr=input('good?',"s");
            if strr=='y'
                preserve_idxs{i,j}(k)=1;
            else
            end
        end
    end
end
