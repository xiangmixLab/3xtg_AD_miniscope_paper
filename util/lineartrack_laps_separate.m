function [nnew,behavnew,laps_dir]=lineartrack_laps_separate(neuron,behav)

behavpos=behav.position;
behavtime=behav.time;

ntime = double(neuron.time);
ntime = ntime(1:2:end);
ntime = resample(ntime,size(neuron.C,2),length(ntime));

max_behavpos_x=abs(max(behavpos(:,1))-min(behavpos(:,1)));
max_behavpos_y=abs(max(behavpos(:,2))-min(behavpos(:,2)));

laps_idx=zeros(size(behavpos,1),1);

if max_behavpos_x>max_behavpos_y % horizontal track
    xpos=behavpos(:,1);
%     xpos=smoothdata(xpos,'SmoothingFactor',0.01);
    % identify laps
    max_lap_iv=max(xpos)-min(xpos);
    max_lap=max(xpos)-max_lap_iv*0.05; % after water removal the end should be clipped
    min_lap=min(xpos)+max_lap_iv*0.05;
    
    xpos=fillmissing(xpos,'linear');
    
    laps_idx(xpos>=max_lap)=1;
    laps_idx(xpos<=min_lap)=2;
    laps_idx1=laps_idx==0;
    
    upos=xpos;
else % vertical track
    ypos=behavpos(:,2);
%     ypos=smoothdata(ypos,'SmoothingFactor',0.01);
    
    % identify laps
    max_lap_iv=max(ypos)-min(ypos);
    max_lap=max(ypos)-max_lap_iv*0.05;
    min_lap=min(ypos)+max_lap_iv*0.05;
    
    ypos=fillmissing(ypos,'linear');
    
    laps_idx(ypos>=max_lap)=1;
    laps_idx(ypos<=min_lap)=2;
    laps_idx1=laps_idx==0;
    
    upos=ypos;
end


laps_list={};
laps_dir=[];

stats=regionprops(laps_idx1);
ctt=1;
for i=1:length(stats)
    if stats(i).Area>30 && floor(stats(i).BoundingBox(2))>1 && floor(stats(i).BoundingBox(2)+stats(i).BoundingBox(4))<size(behavpos,1)% Fs, 1sec, if the 'laps' start from first frame or end last frame, it probably not finished...
        bbox_edge=floor([stats(i).BoundingBox(2),stats(i).BoundingBox(2)+stats(i).BoundingBox(4)]);
        
        if abs(max(upos([bbox_edge(1):bbox_edge(2)]))-min(upos([bbox_edge(1):bbox_edge(2)])))>=max_lap_iv*0.8 % delete incomplete laps
            laps_list{ctt,1}=[bbox_edge(1):bbox_edge(2)];
            if laps_idx(bbox_edge(1)-1)==1 || laps_idx(bbox_edge(1))==1 %right, far side to close side
                laps_dir(ctt,1)=1;
            end
            if laps_idx(bbox_edge(1)-1)==2 || laps_idx(bbox_edge(1))==2%left, close side to far side
                laps_dir(ctt,1)=2;
            end
            ctt=ctt+1;
        end
    end
end

% check
% laps_idx2=laps_idx1*0;
% for i=1:length(laps_list)
% laps_idx2(laps_list{i})=1000;
% end
% plot(laps_idx2);hold on;plot(ypos)

nnew=[];
behavnew=[];
for i=1:length(laps_list)
    idx=zeros(size(behavpos,1),1);
    idx(laps_list{i})=1;
    idx=logical(idx);   
    
    behavnew.position{i,laps_dir(i,1)}=behavpos(idx,:);
    behavnew.time{i,laps_dir(i,1)}=behavtime(idx);
    behavnew.ROI{i,laps_dir(i,1)}=behav.ROI;
    
    idx_resample=interp1(behavtime,double(idx),ntime,'nearest');
    idx_resample1=interp1(behavtime,double(idx),neuron.time,'nearest');

    idx_resample(isnan(idx_resample))=0;
    idx_resample1(isnan(idx_resample1))=0;

    idx_resample=logical(idx_resample);
    idx_resample1=logical(idx_resample1);

    nnew.C{i,laps_dir(i,1)}=neuron.C(:,idx_resample);
    nnew.C_raw{i,laps_dir(i,1)}=neuron.C_raw(:,idx_resample);
    nnew.S{i,laps_dir(i,1)}=neuron.S(:,idx_resample);
    nnew.time{i,laps_dir(i,1)}=neuron.time(idx_resample1);
end

 


