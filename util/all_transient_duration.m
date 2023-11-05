function [medianLength4,maxL4,meanL4,idx4,transientss]=all_transient_duration(fname4,conditions,Fs)

medianLength4=[];
maxL4=[];

idx4=[];
transientss=[];
ctt=1;
for i=1:length(fname4)
    
    load(fname4{i});
    
    % transient length
    for j=1:length(conditions)
        nC=neuronIndividuals_new{conditions(j)}.C;
        
        nC(sum(nC,2)==0,:)=[];
        nC(min(nC,[],2)>5,:)=[];
        nS=neuronIndividuals_new{conditions(j)}.S;

        nStd=std(nS,[],2);
        for k=1:size(nC,1)
            nCk=nC(k,:);
            nCk=nCk-min(nCk(nCk>0));
            nCk(nCk<=nStd(k)*3)=0;
            nC(k,:)=nCk;
        end
        
        nC(sum(nC,2)==0,:)=[];
        nC(min(nC,[],2)>5,:)=[];      
        for k=1:size(nC,1)
            nCk=nC(k,:);
%             [Areas,transientssCur]=findTransients_hwhh(nCk);
            [Areas,transientssCur]=findTransients(nCk);

            if ~isempty(Areas) % abnormal neurons will have no transient area recorded
                medianLength4(ctt,1)=median([Areas])/Fs; % Fs=15
                maxL4(ctt,1)=max([Areas])/Fs; % Fs=15
                meanL4(ctt,1)=mean([Areas])/Fs; % Fs=15
                idx4(ctt,:)=[i,conditions(j)]; % mouse, condition
                transientss=[transientss;transientssCur];
                ctt=ctt+1;
            end
        end
    end
end        

