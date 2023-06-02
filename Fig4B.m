rng('shuffle');
tmax = 5000;

% parameters
ntot = 2000;
rho = 1e10; % (aa/um^3)
a = 20;
g = ones(ntot,1);
g(2) = 5; % g_r
L = 500*ones(ntot,1);
L(1) = 1e3; % L_n
L(2) = 1e4; % L_r
K_nmean = 6000; % <K_{n,i}> (um^(-3)), 6000 um^(-3)~10 umol/L
CVK_n = 0; % CV of K_{n,i}
sigma = sqrt(log(1+CVK_n^2)); % for X~LogNormal(mu,sigma^2), <X>=e^(mu+sigma^2/2), CV=\sqrt{e^(sigma^2)-1}
mu = log(K_nmean)-sigma^2/2;
K_n = exp(mu+sigma*randn(ntot,1)); % sampling K_{n,i}~LogNormal(mu,sigma^2)
K_n(1) = K_nmean;
K_n(2) = K_nmean;
gx = 10; % copy number of activator and inhibitor
K_n(3:2+gx) = 12000;
K_n(3+gx:2+2*gx) = 1000;

v_n = 720; % (aa/min)
tau_m = 5*ones(ntot,1); % (min)
Gamma_r = 10*ones(ntot,1); % (min^(-1))
K_r = 6000*ones(ntot,1); % (um^(-3))
v_r = 720; % (aa/min)
tau_p = 1e4*ones(ntot,1); % (min)
tau_p(3:2+2*gx) = 5;

T_D = 120; % doubling time (min)
mu = log(2)/T_D;
n_c = 1e4;
Gamma_n = zeros(ntot,1);
Gamma_n(2) = mu*v_n*(n_c-sum(g))*(1+Gamma_r(3)*L(3)/v_r)/(g(2)*(Gamma_r(3)*L(3)-mu*(1.1*L(3)-0.1*L(1)-L(2)))); % \Gamma_{n,r}
Gamma_n(1) = 0.1*Gamma_n(2)*g(2)/g(1); % \Gamma_{n,n}
Gamma_n(3:end) = (v_n*(n_c-sum(g))-Gamma_n(1)*g(1)*L(1)-Gamma_n(2)*g(2)*L(2))/(L(3)*sum(g(3:end))); % \Gamma_{n,i}, i=3,...,ntot

% initial conditions
V_0 = 0.1;
c_n = a*rho*g(1)*Gamma_n(1)/sum(Gamma_n.*g.*L);
n_0 = c_n*V_0/a;
p_0 = Gamma_n.*g*n_0/(Gamma_n(1)*g(1));
F_n = fzero(@(F_n) n_0*(1-F_n)-sum(g.*(1+Gamma_n.*L./v_n)*F_n*c_n./(F_n*c_n+K_n)),[0 1]);
m_0 = Gamma_n.*g*F_n*c_n./(F_n*c_n+K_n).*tau_m; % initial condition m_i=k_{n,i}tau_{m,i}

dd = 0.1;
genmax = 5000; % the number of generations
threshold = 0.45;

% WT
[meanVb0,CVVb0,Vb0] = stochastictrace1(threshold,dd,genmax,m_0,p_0,tmax,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p,gx);

% delete one activator
g(3) = 0;
m_0(3) = 0;
p_0(3) = 0;
[meanVb,CVVb,Vb] = stochastictrace1(threshold,dd,genmax,m_0,p_0,tmax,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p,gx);