function ratemap_plot(fr,countTime,smoothsign,colorbar_sign,caxis_range)
if smoothsign==1
    fr=filter2DMatrices(fr,1);
end
fr(countTime==0)=nan;
fr_pad=zeros(size(fr,1)+1,size(fr,2)+1);
fr_pad(1:size(fr,1),1:size(fr,2))=fr;
pcolor(fr_pad); % pcolor seems ignore the last row and col

if colorbar_sign
    colorbar;
end

if ~isempty(caxis_range)
    caxis(caxis_range);
end
colormap(jet);
axis ij
shading flat;
axis image
axis off