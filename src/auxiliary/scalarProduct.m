function Unorm=scalarProduct(Uface,fNormal);

Unorm=sum(Uface.data .* fNormal,1);

end