% if no fs2, please left it empty
function rateMap_corr=rateMap_correlation(fs1,fs2,ct1,ct2,resize_sign,smooth_sign)
rateMap_corr=[-1];
ctt=1;

if ~isempty(fs2)
    for i=1:length(fs1)
        if ~isempty(fs1{i})&&~isempty(fs2{i})
            cell1=fs1{i};
            cell2=fs2{i};
            if smooth_sign==1
                cell1=filter2DMatrices(cell1,1);
                cell2=filter2DMatrices(cell2,1);
            end
            if resize_sign==1
                cell2=imresize(cell2,[size(cell1,1),size(cell1,2)]);
                ct1=imresize(ct1,[size(cell1,1),size(cell1,2)]);
                ct2=imresize(ct2,[size(ct1,1),size(ct1,2)]);
            end
%             cell1(ct1==0)=nan; % get rid of non trespassing periods, avoid unnecessary correlation boost
%             cell2(ct2==0)=nan;
            cell1t=reshape(cell1,size(cell1,1)*size(cell1,2),1);
            cell2t=reshape(cell2,size(cell2,1)*size(cell2,2),1);
    %         idx=~isnan(cell1t);
            if isempty(cell1t)
                cell1t=cell2t*0;
            end
            if isempty(cell2t)
                cell2t=cell1t*0;
            end
            cell1t(isnan(cell1t))=0;
            cell2t(isnan(cell2t))=0;
            
%             idx0=unique([find(cell1t==0);find(cell2t==0)]);
%             cell1t(idx0)=[];
%             cell2t(idx0)=[]; % remoce 0 to reduce unnecessary corr
% 
            rm_corr=corrcoef(cell1t,cell2t);
            try
                rateMap_corr(ctt)=rm_corr(1,2);
            catch
                rateMap_corr(ctt)=nan;
            end
            ctt=ctt+1;
        else
            rateMap_corr(ctt)=0;
            ctt=ctt+1;
        end
     end
else
    if ~isempty(fs1)
        for i=1:length(fs1)-1
            for j=i+1:length(fs1)
                cell1=fs1{i};
                cell2=fs1{j};
                cell1=filter2DMatrices(cell1,1);
                cell2=filter2DMatrices(cell2,1);
                cell1(ct1==0)=nan; % get rid of non trespassing periods, avoid unnecessary correlation boost
                cell2(ct1==0)=nan;
                cell1=reshape(cell1,size(cell1,1)*size(cell1,2),1);
                cell2=reshape(cell2,size(cell2,1)*size(cell2,2),1);
                if isempty(cell1)
                    cell1=cell2*0;
                end
                if isempty(cell2)
                    cell2=cell1*0;
                end
                cell1(isnan(cell1))=0;
                cell2(isnan(cell2))=0;

                idx0=unique([find(cell1==0);find(cell2==0)]);
                cell1(idx0)=[];
                cell2(idx0)=[]; % remoce 0 to reduce unnecessary corr

                rm_corr=corrcoef(cell1,cell2);
                rateMap_corr(ctt)=rm_corr(1,2);
                ctt=ctt+1;
            end
        end
    else
        rateMap_corr=[-1];
    end
end