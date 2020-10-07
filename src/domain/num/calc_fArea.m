function fArea = calc_fArea(nF,fNbVLoc,fNbV,vCoord)
%for 2D only
    fArea=zeros(1,nF);
    
    for i=1:nF
      Vind=fNbV(fNbVLoc*(i-1)+1:fNbVLoc*i); % indicies of neighbour verticies
      fArea(i)=pdist([vCoord(:,Vind(1))';vCoord(:,Vind(2))'],'euclidean');   
    end
end



