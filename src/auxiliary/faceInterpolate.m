function Uface=faceInterpolate(dom,U)
tempUface=zeros(2,dom.nF);
for i=1:dom.nF
  PC=dom.fNbC(2*i-1); %physical cell index
  NBC=dom.fNbC(2*i); %neighbor cell
  tempUface(:,i)=(U.data(:,PC)+U.data(:,NBC))/2.0;
end

Uface = Field(dom.allFaces,1);     % Velocity [m/s] (vector);
set(Uface,tempUface);

end