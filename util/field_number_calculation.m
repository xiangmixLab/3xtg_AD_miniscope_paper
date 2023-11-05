function numFields=field_number_calculation(firingrate)

numFields=[];
for i=1:length(firingrate)
    ft=firingrate{i};
    
    if ~isempty(ft)
        ft=filter2DMatrices(ft,1);
        ft=imadjust(adapthisteq(ft));
        ft(ft<max(ft(:))*0.1)=0;
        ft=ft.*ft;
        ft=ft>max(ft(:))*0.2;
%         ft=bwareaopen(ft,9);
        
        stats=regionprops(ft);
        
        numFields(i)=length(stats);
    else
        numFields(i)=0;
    end
end