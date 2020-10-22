function fieldFace=faceInterpolate(dom,field)
tempFace=zeros(field.dim,dom.nF);
for i=1:dom.nF
  PC=dom.fNbC(2*i-1); %physical cell index
  NBC=dom.fNbC(2*i); %neighbor cell
  fXiLambda=dom.fXiLambda(i);
  tempFace(:,i)=( (1-fXiLambda)*field.data(:,PC) + fXiLambda*field.data(:,NBC) );
end

fieldFace = Field(dom.allFaces,1);     % Velocity [m/s] (vector);
set(fieldFace,tempFace);

end