function [nnew1,bnew1]=rearrange_laps(nnew,behavnew)

nnew1={};
bnew1={};
for i=1:size(nnew,1)
    for j=1:size(nnew,2)
        nC=nnew{i,j}.C;
        nS=nnew{i,j}.S;
        nCraw=nnew{i,j}.C_raw;
        nTime=nnew{i,j}.time;
        
        nnew_tmp={};
        for k1=1:size(nC,1)
            for k2=1:size(nC,2)
                if ~isempty(nC{k1,k2})
                    nnew_tmp{k1,k2}.C=nC{k1,k2};
                    nnew_tmp{k1,k2}.S=nS{k1,k2};
                    nnew_tmp{k1,k2}.C_Raw=nCraw{k1,k2};
                    nnew_tmp{k1,k2}.time=nTime{k1,k2};
                end
            end
        end
        
        
        position=behavnew{i,j}.position;
        btime=behavnew{i,j}.time;
        behav_tmp={};
        for k1=1:size(position,1)
            for k2=1:size(position,2)
                if ~isempty(position{k1,k2})
                    behav_tmp{k1,k2}.position=position{k1,k2};
                    behav_tmp{k1,k2}.time=btime{k1,k2};
                    behav_tmp{k1,k2}.ROI=behavnew{i,j}.ROI{k1,k2};
                end
            end
        end
        
        nnew1{i,j}=nnew_tmp;
        bnew1{i,j}=behav_tmp;
    end
end