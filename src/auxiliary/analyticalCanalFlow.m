function u = analyticalCanalFlow(ymesh,dp,H,mu,Utop)

u = dp*H^2/(2*mu)*( ymesh.^2/H^2 - ymesh/H ) + Utop * ymesh/H;

end