function fXi = calc_fXi(nF,cCoord,fNbCLoc,fNbC)
%2D only
  fXi=zeros(2,nF);
    
    for i=1:nF
      Cind=fNbC(fNbCLoc*(i-1)+1:fNbCLoc*i); % indicies of neighbour cells
      fXi(:,i)= cCoord(:,Cind(2))-cCoord(:,Cind(1));
    end

end



