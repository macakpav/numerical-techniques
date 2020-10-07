function fNormal = calc_fNormal(fTangent)
%for 2D only
    fNormal=zeros(size(fTangent));
    fNormal(1,:)=-fTangent(2,:); 
    fNormal(2,:)=fTangent(1,:);
end