function field = analyticalHotPlateDiffusion(dom, T1, T2, H, W)
field = zeros(1,dom.nC);
n_max=800;

for i=1:dom.nPc
  x=dom.cCoord(1,i);
  y=dom.cCoord(2,i);
  sum=0;
  
  for n = n_max:-1:1
    k=((-1)^(n+1)+1)/n;
    a = sin(n*pi*x/W);
    b = sinh(n*pi*y/W);
    c = sinh(n*pi*H/W);
    
    add = k*a*b/c;
    
    if add ~= add % is NaN
      if c == Inf
        if a == Inf
          add= k * exp(n*pi*x/W+n*pi*y/W-n*pi*H/W);
        else
          add= k * a * exp(n*pi*y/W-n*pi*H/W);
        end
      else
      end
    end
    sum=sum+add;
  end
  field(i)=T1+(T2-T1)*2/pi*sum;
  
end




end