function neuronIndividuals=neuronIndividuals_new_split(neuron,frames,interval)

startt=1;
startt1=1;
neuronIndividuals={};
for i=1:floor((size(neuron.C,2)-frames)/interval)+1
    neuronIndividuals{i}.C=neuron.C(:,startt:min(startt+interval,size(neuron.C,2)));
    neuronIndividuals{i}.S=neuron.S(:,startt:min(startt+interval,size(neuron.C,2)));
    neuronIndividuals{i}.C_df=[];
    neuronIndividuals{i}.trace=[];
    neuronIndividuals{i}.C_raw=neuron.C_raw(:,startt:min(startt+interval,size(neuron.C,2)));
    neuronIndividuals{i}.time=neuron.time(startt1:min(startt1+interval*2,length(neuron.time)));
    
    disp(num2str([startt,startt+interval]))
    startt=startt+round(interval);
    startt1=startt1+round(2*interval);
end
