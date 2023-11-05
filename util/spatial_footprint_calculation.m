function [footprint,centroid]=spatial_footprint_calculation(neuron,varargin)

if ~isfield(neuron,'imageSize')||isempty(neuron.imageSize)
    neuron.imageSize=size(neuron.Cn);
end


d1=neuron.imageSize(1);
d2=neuron.imageSize(2);


summ=zeros(d1*d2,1);

centroid=[];

if length(varargin)>0
    thresh=varargin{1};
else
    thresh=0;
end

for i=1:size(neuron.A,2)
    Ai=full(neuron.A(:,i));

    Ai=Ai*1/max(Ai(:));
    ind=find(Ai>thresh);%the pixels at which a neuron shows up
    summ(ind)=Ai(ind);
    Ait=reshape(Ai,d1,d2);
    Ait=Ait>thresh;
    stats=regionprops(Ait,'Centroid');
    centroid(i,:)=stats.Centroid;
%     centroid(i,2)=d1-centroid(i,2); % ij to xy
end

footprint=reshape(summ,d1,d2);

ft_bin=footprint>0;

