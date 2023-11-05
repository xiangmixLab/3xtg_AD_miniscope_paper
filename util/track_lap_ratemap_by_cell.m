function [lap_ratemap,lap_ratemap_full]=track_lap_ratemap_by_cell(fr_all,cellNum)

lap_ratemap={};
lap_ratemap_full={};

sz1=[];
sz2=[];
for i=1:size(fr_all,1)
    % find the first non-empty to determine size
    if ~isempty(fr_all{i,1})
        for j=1:length(fr_all{i,1})
            if ~isempty(fr_all{i,1}{j})
                sz1=size(fr_all{i,1}{j});
            end
        end
    end
    if ~isempty(fr_all{i,2})
        for j=1:length(fr_all{i,2})
            if ~isempty(fr_all{i,2}{j})
                sz2=size(fr_all{i,2}{j});
            end
        end
    end
    if ~isempty(sz1)&&~isempty(sz2)
        break;
    end
end

if sz1(1)<sz1(2) % hor condition
    lap_lmap=[];
    lap_rmap=[];
    
    ctt1=1;
    ctt2=1;
    
    for i=1:size(fr_all,1)    
        if ~isempty(fr_all{i,1})&&~isempty(fr_all{i,1}{cellNum})
            lap_lmap(ctt1,:)=mean(filter2DMatrices(fr_all{i,1}{cellNum},1),1);
            ctt1=ctt1+1;
        end
        if ~isempty(fr_all{i,2})&&~isempty(fr_all{i,2}{cellNum})
            lap_rmap(ctt2,:)=mean(filter2DMatrices(fr_all{i,2}{cellNum},1),1);
            ctt2=ctt2+1;
        end
    end

    lap_ratemap{1}=lap_lmap;
    lap_ratemap{2}=lap_rmap;
else % vert condition
    lap_lmap=[];
    lap_rmap=[];

    ctt1=1;
    ctt2=1;
    
    for i=1:size(fr_all,1)    
        if ~isempty(fr_all{i,1})&&~isempty(fr_all{i,1}{cellNum})
            lap_lmap(ctt1,:)=mean(filter2DMatrices(fr_all{i,1}{cellNum},1),2)';
            ctt1=ctt1+1;
        end
        if ~isempty(fr_all{i,2})&&~isempty(fr_all{i,2}{cellNum})
            lap_rmap(ctt2,:)=mean(filter2DMatrices(fr_all{i,2}{cellNum},1),2)';
            ctt2=ctt2+1;
        end
    end

    lap_ratemap{1}=lap_lmap;
    lap_ratemap{2}=lap_rmap;
end