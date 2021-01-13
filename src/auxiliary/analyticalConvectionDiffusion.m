function phi = analyticalConvectionDiffusion(x_coords, L, u, rho, k, phi_0, phi_L)

phi = zeros(size(x_coords));
i = 1;

if u == 0 % pure diffusion
  for x = x_coords
    phi(i) = phi_0 + (phi_L-phi_0) / L * x;
    i = i + 1;
  end

else % convection-diffusion
  P = u * L * rho / k;
  for x = x_coords
    phi(i) = phi_0 + (phi_L - phi_0) * (exp(P * x / L) - 1) / (exp(P) - 1);
    i = i + 1;
  end
end

end