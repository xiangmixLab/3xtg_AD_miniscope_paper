function [firingRateAll,countAll,countTimeAll,countTime,amplitudeAll,amplitude_rateAll,binInfo,countTime_pt] = calculatingCellSpatialForSingleData_040321(neuron,behavpos,behavtime,behavROI,binsize,segments,threshold,temp,intensity,plotting,countTimeThresh,small_velo)
%% Inputs:
%     neuron: a source2D variable, including identified neurons with traces and spatial information, which is obtained by runing cnmfe codes
%     behav: behavior information, which is obtained by using Tristan's code
% Segment: a vector, e.g,1:10 (display the traces of the first 10 identified neurons)
% threshold:the threshold above which the neuron is active, e.g.,0.1
% temp:judge whether the neuron is active using neuron.trace or neuron.S; temp = 'trace' ot temp = 'S'
% downsampling: if downsampling = true, then do downsampling for neuron.time; default is false
% intensity: if intensity = true, then diplay the peak map; otherwise display the heatmap of firating rate; default is false
%%% example usage: e.g.1. plottingCellSpatialForSingleData(neuron,behav,1:5)
%%%                e.g.2, plottingCellSpatialForSingleData(neuron,behav,1:5,0.1,'trace',true,true)
%%%                e.g.3, plottingCellSpatialForSingleData(neuron,behav,1:5,0.1,'trace')
if ~exist('plotting','var') || isempty(plotting)
    plotting = false;
end
if ~exist('intensity','var') || isempty(intensity)
    intensity = false;
end
if ~exist('temp','var') || isempty(temp)
    temp = 'S';
end
if ~exist('threshold','var') || isempty(threshold)
    threshold = 0.1;
end
if ~exist('segments','var') || isempty(segments)
    segments = 1:size(neuron.C,1);
end
if ~exist('binsize','var') || isempty(binsize)
    binsize = 15;%attention: value has been 5, 10, 15
end
if ~exist('countTimeThresh','var') || isempty(countTimeThresh)
    countTimeThresh = [0 inf]; % unit: sec
end

if ~exist('small_velo','var') || isempty(small_velo)
    small_velo = -1; % unit: mm/sec
end
% behavpos = behav.position;
% behavtime = behav.time;
% behavROI = behav.ROI;

downsampling = length(neuron.time)/size(neuron.C,2);
if downsampling ~= 1
    %     downsampling == 2
    neuron.time = double(neuron.time);
    neuron.time = neuron.time(1:downsampling:end);
    neuron.time = resample(neuron.time,size(neuron.C,2),length(neuron.time));
end
t = find(diff(behavtime)<=0);
while ~isempty(t)
    behavtime(t+1) = behavtime(t)+1;
    t = find(diff(behavtime)<=0);
end
behavpos = fillmissing(behavpos,'nearest');
neuron.pos = interp1(behavtime,behavpos,neuron.time); %%
neuron.pos = fillmissing(neuron.pos,'nearest');
folderName = 'FiguresCellSpatial';
% if ~exist(folderName,'dir')
%     mkdir(folderName)
% end

fpath=[folderName];
global ts pos1 pos2
num = length(segments);

pos1 = 0:binsize:ceil(behavROI(3));
pos2 = 0:binsize:ceil(behavROI(4));

% pos1 = 0:binsize:ceil(max(neuron.pos(:,1)));
% pos2 = 0:binsize:ceil(max(neuron.pos(:,2)));
% pos1 = 0:binsize:round(max(neuron.pos(:,1)));
% pos2 = 0:binsize:round(max(neuron.pos(:,2)));
% pos1 = 0:binsize:round(max(neuron.pos(:,1)));
% pos2 = 0:binsize:round(max(neuron.pos(:,2)));
% pos1 = floor(min(neuron.pos(:,1))):binsize:ceil(max(neuron.pos(:,1)));
% pos2 = floor(min(neuron.pos(:,2))):binsize:ceil(max(neuron.pos(:,2)));
% decrement1 = max(neuron.pos(:,1)-min(neuron.pos(:,1)))*0.005;
% decrement2 = max(neuron.pos(:,2)-min(neuron.pos(:,2)))*0.005;
% pos1 = floor(min(neuron.pos(:,1))+decrement1):binsize:ceil(max(neuron.pos(:,1))-decrement1);
% pos2 = floor(min(neuron.pos(:,2))+decrement2):binsize:ceil(max(neuron.pos(:,2))-decrement2);

% sprintf("The number of bins is %dx%d",length(pos1),length(pos2))

% if max(pos1) < floor(max(neuron.pos(:,1)))
%     pos1 = [pos1 max(pos1)+binsize];
% end
% if max(pos2) < floor(max(neuron.pos(:,2)))
%     pos2 = [pos2 max(pos2)+binsize];
% end
binInfo.binsize = binsize;binInfo.pos1 = pos1;binInfo.pos2 = pos2;
binInfo.xpos = pos1;
binInfo.ypos = pos2;
% pos1 = 0:binsize:ceil(behavROI(:,3));%we changed to ROI
% pos2 = 0:binsize:ceil(behavROI(:,4));

%% small velo determine % added 061219
if small_velo>0
    d_behavpos=zeros(size(behavpos,1)-1,1);
    for ll=2:size(behavpos,1)
        d_behavpos(ll-1,1)=norm((behavpos(ll,:)-behavpos(ll-1,:)));
    end
    d_behavtime=diff(behavtime)/1000;
    velo=d_behavpos./d_behavtime;
    small_velo_idx=find(velo<small_velo);
    slow_period=behavtime*0;
    slow_period(small_velo_idx+1)=1;
end
%     countTime_slow_period=calculateSlowPeriodcountTime(behavpos,behavtime,slow_period,binsize,behavROI);
%     countTime_slow_period_idx=countTime_slow_period>0;
%     
%     countTime=countTime.*(~countTime_slow_period_idx);


%% countTime calculation
countTime = zeros(length(pos1),length(pos2));
countTime_pt = zeros(length(pos1),length(pos2));
ts = 1;
behavpos_bin=ceil(behavpos./binsize);
behavpos_bin(behavpos_bin<=0)=1;

behavpos_bin(behavpos_bin(:,1)>length(pos1),1)=length(pos1);
behavpos_bin(behavpos_bin(:,2)>length(pos2),2)=length(pos2);

diff_behavtime=diff(behavtime);
for i=1:size(diff_behavtime,1)
    countTime(behavpos_bin(i,1),behavpos_bin(i,2))=countTime(behavpos_bin(i,1),behavpos_bin(i,2))+diff_behavtime(i);
    countTime_pt(behavpos_bin(i,1),behavpos_bin(i,2))=countTime_pt(behavpos_bin(i,1),behavpos_bin(i,2))+1;
end
countTime = countTime'/1000;%purpose: because behavtime is recorded in milisec.
countTime_pt = countTime_pt';


%%
%maxCount = zeros(1,num);
maxCount = 5; %add this line when you obtain the maximum of maxCount after running once
% firingRateAll = zeros(length(pos1),length(pos2),num);
% countAll = zeros(length(pos1),length(pos2),num);
firingRateAll = cell(1,num);
countAll = cell(1,num);
countTimeAll = cell(1,num);
amplitudeAll = cell(1,num);
amplitude_rateAll=cell(1,num);
for k = 1:num
    if strcmpi(temp,'trace')
        if length(threshold) <= 1&&threshold==0.1
            %         thresh = (max(neuron.trace(Segment(k),:))-min(neuron.trace(Segment(k),:)))*threshold; % the threshold above which the neuron is active
            thresh = (max(neuron.trace(segments(k),:))-0)*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(segments(k));
        end
        %         idx = find(neuron.trace(segments(k),:)>thresh);
        [pks,locs] = findpeaks(neuron.trace(segments(k),:),'MinPeakHeight',thresh);
        idx = locs;
    elseif strcmpi(temp,'S')
        if length(threshold) <= 1&&threshold==0.1
            %         thresh = (max(neuron.S(Segment(k),:))-min(neuron.S(Segment(k),:)))*threshold; % the threshold above which the neuron is active
            thresh = (max(neuron.S(segments(k),:))-0)*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(segments(k));
        end
        idx = find(neuron.S(segments(k),:)>thresh);
    elseif strcmpi(temp,'C')
        if length(threshold) <= 1&&threshold==0.1
            %         thresh = (max(neuron.C(Segment(k),:))-min(neuron.C(Segment(k),:)))*threshold; % the threshold above which the neuron is active
            thresh = (max(neuron.C(segments(k),:))-0)*threshold; % the threshold above which the neuron is active
        else
            thresh = threshold(segments(k));
        end
        idx = find(neuron.C(segments(k),:)>thresh);
        %         [pks,locs,w,p] = findpeaks(neuron.C(segments(k),:),'MinPeakDistance',100,'MinPeakHeight',thresh);
%         [pks,locs] = findpeaks(neuron.C(segments(k),:),'MinPeakHeight',thresh);
%         idx = locs;
    end
    if ~isempty(idx)
        
        if small_velo>0 % added 061219
            slow_periodt=interp1(behavtime,slow_period,neuron.time,'next'); %%resample(slow_period,size(neuron.C,2),length(slow_period));
%             slow_periodt=resample(slow_period,size(neuron.C,2),length(slow_period)); %%resample(slow_period,size(neuron.C,2),length(slow_period));
            idx(ismember(idx,intersect(idx,find(slow_periodt==1))))=[];
        end
        
        
        count = countingFiringBins(idx,neuron,binsize);
        amplitude = countingAmplitudeFiringBins(idx,neuron,segments(k),binsize);
        %% small time threshold && long time threshold
%         countTime_smaller_than_thr=countTime<countTimeThresh(1); %0.2sec
%         count(countTime_smaller_than_thr)=0;
%         amplitude(countTime_smaller_than_thr) = 0;
%         countTime_larger_than_thr=countTime>countTimeThresh(2); %10sec
%         count(countTime_larger_than_thr)=0;
%         amplitude(countTime_larger_than_thr) = 0;
        countTime(countTime<countTimeThresh(1))=0;
        countTime(countTime>countTimeThresh(2))=0;
        %         amplitude = amplitude./count; % change amplitude here
        count(countTime == 0) = 0;
        amplitude(countTime == 0) = 0;
        amplitude_rate = amplitude./countTime;
        amplitude_rate(countTime == 0) = 0;
        
        firingRate = count./countTime;
        firingRate(countTime == 0) = 0;
        %  firingRate(isnan(firingRate))=0;
        
        firingRate2 = nan(size(firingRate)+1);
        firingRate2(1:end-1,1:end-1) = firingRate;
        countTime2 = nan(size(countTime)+1);
        countTime2(1:end-1,1:end-1) = countTime;
        firingRate2(countTime2 == 0) = 0;
        %  firingRate2(isnan(firingRate2))=0;
        %         maxCount = max(firingRate(:));
        if plotting
            figure;
            %              pcolor(firingRate2);
            firingRateSmoothing = filter2DMatrices(firingRate2, 1);
            pcolor(firingRateSmoothing);
            colormap(jet)
            %  maxCount(k) = max(firingRate(:)); % comment this line when you obtain the maximum of maxCount after running once
            caxis([0,max(maxCount)])
            colorbar;
            set(gca, 'color', 'w', 'ydir', 'reverse')
            shading flat;
            if intensity
                shading interp
            end
            axis image
            axis off;
            axis ij
            title(['Cell #', num2str(segments(k))],'FontName','Arial','FontSize',10,'FontWeight','bold')
            hold off
            saveas(gcf,fullfile(fpath,['CellSpatialMatch',num2str(segments(k)),'.tif']))
            saveas(gcf,fullfile(fpath,['CellSpatialMatch',num2str(segments(k)),'.fig']))
        end
        
        firingRateAll{k} = firingRate;
        countAll{k} = count;
        countTimeAll{k} = countTime;
        amplitudeAll{k} = amplitude;
        amplitude_rateAll{k} = amplitude_rate;
        %     firingRateAll(:,:,k) = firingRate;
        %     countAll(:,:,k) = count;
    end
    
    %max(maxCount)
    
end

function count=countingFiringBins(idx,neuron,binsize)
global pos1 pos2
count = zeros(length(pos1),length(pos2));
for i = 1:length(idx)
    % if idx(i)<=length(neuron.pos) %sometimes in behav data, the position is less longer than neuron.S, which cause crash
%     [~,idxx] = find(pos1 <= neuron.pos(idx(i),1), 1, 'last');
%     [~,idyy] = find(pos2 <= neuron.pos(idx(i),2), 1, 'last');
    
%     count(idxx,idyy) = count(idxx,idyy)+1;
      coorr=[ceil(neuron.pos(idx(i),1)/binsize),ceil(neuron.pos(idx(i),2)/binsize)];
      coorr(coorr<=0)=1; % left boundary readjust
      coorr(1)=min(length(pos1),coorr(1));% right boundary readjust
      coorr(2)=min(length(pos2),coorr(2));% right boundary readjust
      
      count(coorr(1),coorr(2)) = count(coorr(1),coorr(2))+1;
    % end
end
count = count';
% count(end+1,:) = count(end,:);
% count(:,end+1) = count(:,end);

function amplitude = countingAmplitudeFiringBins(idx,neuron,k,binsize)
global pos1 pos2
amplitude = zeros(length(pos1),length(pos2));
for i = 1:length(idx)
    % if idx(i)<=length(neuron.pos) %sometimes in behav data, the position is less longer than neuron.S, which cause crash
%     [~,idxx] = find(pos1 <= neuron.pos(idx(i),1), 1, 'last');
%     [~,idyy] = find(pos2 <= neuron.pos(idx(i),2), 1, 'last');
%     amplitude(idxx,idyy) = amplitude(idxx,idyy)+neuron.C(k,idx(i));
     coorr=[ceil(neuron.pos(idx(i),1)/binsize),ceil(neuron.pos(idx(i),2)/binsize)];
     coorr(coorr<=0)=1;
     coorr(1)=min(length(pos1),coorr(1));
     coorr(2)=min(length(pos2),coorr(2));
      

     amplitude(coorr(1),coorr(2)) = amplitude(coorr(1),coorr(2))+neuron.C(k,idx(i));
    % end
end
amplitude = amplitude';

