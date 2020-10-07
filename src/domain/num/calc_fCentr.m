function fCentr = calc_fCentr(nF,fNbVLoc,fNbV,vCoord)
    fCentr=zeros(2,nF);
    
    for i=1:nF
      Vind=fNbV(2*i-1:2*i); % indicies of neighbour verticies
      fCentr(:,i)=(vCoord(:,Vind(1))+vCoord(:,Vind(2)))/2;
    end

end



