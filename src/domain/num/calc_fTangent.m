function fTangent = calc_fTangent(fNbVLoc,fNbV,vCoord,fArea)
%for 2D only
    fTangent=zeros(2,length(fArea));
    
    for i=1:length(fArea)
      Vind=fNbV(fNbVLoc*(i-1)+1:fNbVLoc*i); % indicies of neighbour verticies
      fTangent(:,i)= (vCoord(:,Vind(2))-vCoord(:,Vind(1)))* (1./fArea(i));
    end
    
end



