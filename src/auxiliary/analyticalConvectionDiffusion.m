function phi=analyticalConvectionDiffusion(x_coords, L, u, rho, k, phi_0, phi_L)
  phi=zeros(size(x_coords));
  P=u*L*rho/k;
  i=1;
  for x = x_coords
    phi(i)=phi_0+(phi_L-phi_0)*(exp(P*x/L)-1)/(exp(P)-1);
    i=i+1;
  end
end