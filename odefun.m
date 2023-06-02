function dY = odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p)
dY = zeros(2*ntot,1);
M = sum(Y(ntot+1:2*ntot).*L);
c_n = rho*a*Y(ntot+1)/M;
F_n = fzero(@(F_n) Y(ntot+1)*(1-F_n)-sum(g.*(1+Gamma_n.*L./v_n)*F_n*c_n./(F_n*c_n+K_n)),[0 1]);
dY(1:ntot) = Gamma_n.*g*F_n*c_n./(F_n*c_n+K_n)-Y(1:ntot)./tau_m;
c_r = rho*Y(ntot+2)/((1-1/a)*M);
F_r = fzero(@(F_r) Y(ntot+2)*(1-F_r)-sum(Y(1:ntot).*(1+Gamma_r.*L./v_r)*F_r*c_r./(F_r*c_r+K_r)),[0 1]);
dY(ntot+1:2*ntot) = Gamma_r.*Y(1:ntot)*F_r*c_r./(F_r*c_r+K_r)-Y(ntot+1:2*ntot)./tau_p;
end