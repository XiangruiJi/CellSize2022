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
K_n = 6000*ones(ntot,1); % (um^(-3))
K_n(3) = 12000;
K_n(4) = 1000;
v_n = 720; % (aa/min)
tau_m = 5*ones(ntot,1); % (min)
Gamma_r = 10*ones(ntot,1); % (min^(-1))
K_r = 6000*ones(ntot,1); % (um^(-3))
v_r = 720; % (aa/min)
tau_p = 1e4*ones(ntot,1); % (min)
tau_p(3) = 5;

T_D = 120; % doubling time (min)
mu = log(2)/T_D;
n_c = 1e4;
Gamma_n = zeros(ntot,1);
Gamma_n(2) = mu*v_n*(n_c-sum(g))*(1+Gamma_r(3)*L(3)/v_r)/(g(2)*(Gamma_r(3)*L(3)-mu*(1.1*L(3)-0.1*L(1)-L(2)))); % \Gamma_{n,r}
Gamma_n(1) = 0.1*Gamma_n(2)*g(2)/g(1); % \Gamma_{n,n}
Gamma_n(3:end) = (v_n*(n_c-sum(g))-Gamma_n(1)*g(1)*L(1)-Gamma_n(2)*g(2)*L(2))/(L(3)*sum(g(3:end))); % \Gamma_{n,i}, i=3,...,ntot
a0 = 10;
Gamma_n(4) = Gamma_n(4)/a0;

% initial conditions
V_0 = 0.1;
c_n = a*rho*g(1)*Gamma_n(1)/sum(Gamma_n.*g.*L);
n_0 = c_n*V_0/a;
p_0 = Gamma_n.*g*n_0/(Gamma_n(1)*g(1));
F_n = fzero(@(F_n) n_0*(1-F_n)-sum(g.*(1+Gamma_n.*L./v_n)*F_n*c_n./(F_n*c_n+K_n)),[0 1]);
m_0 = Gamma_n.*g*F_n*c_n./(F_n*c_n+K_n).*tau_m; % initial condition m_i=k_{n,i}tau_{m,i}
Y0 = [m_0;p_0];

threshold = 0.1;
gamma = 0.4;
options = odeset('Events',@(t,Y) myEventsFcn_3(t,Y,ntot,threshold));
t_all = [];
m_super = [];
m_sub = [];
p_super = [];
p_sub = [];
V_all = [];
F_nall = [];
F_rall = [];
tlast = 0;
tend = zeros(8,1);
T_M = 60; % (min)
for generation = 1:8
    [t,Y,~,~,~] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,tmax],Y0,options);
    t_all = [t_all;t+tlast];
    m_super = [m_super;Y(:,3)];
    m_sub = [m_sub;Y(:,4)];
    p_super = [p_super;Y(:,ntot+3)];
    p_sub = [p_sub;Y(:,ntot+4)];
    V = sum(Y(:,ntot+1:2*ntot).*L',2)/rho;
    V_all = [V_all;V];
    for k=1:length(t)
        M = sum(Y(k,ntot+1:2*ntot)'.*L);
        c_n = rho*a*Y(k,ntot+1)/M;
        F_n = fzero(@(F_n) Y(k,ntot+1)*(1-F_n)-sum(g.*(1+Gamma_n.*L./v_n)*F_n*c_n./(F_n*c_n+K_n)),[0 1]);
        F_nall = [F_nall;F_n];
        c_r = rho*Y(k,ntot+2)/((1-1/a)*M);
        F_r = fzero(@(F_r) Y(k,ntot+2)*(1-F_r)-sum(Y(k,1:ntot)'.*(1+Gamma_r.*L./v_r)*F_r*c_r./(F_r*c_r+K_r)),[0 1]);
        F_rall = [F_rall;F_r];
    end
    if t(end)>=tmax
        break;
    end
    Y0 = Y(end,:)';
    tlast = t_all(end);
    % DNA replication
    g = g*2;
    % divide after T_M
    [t,Y] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,T_M],Y0);
    t_all = [t_all;t+tlast];
    m_super = [m_super;Y(:,3)];
    m_sub = [m_sub;Y(:,4)];
    p_super = [p_super;Y(:,ntot+3)];
    p_sub = [p_sub;Y(:,ntot+4)];
    V = sum(Y(:,ntot+1:2*ntot).*L',2)/rho;
    V_all = [V_all;V];
    for k=1:length(t)
        M = sum(Y(k,ntot+1:2*ntot)'.*L);
        c_n = rho*a*Y(k,ntot+1)/M;
        F_n = fzero(@(F_n) Y(k,ntot+1)*(1-F_n)-sum(g.*(1+Gamma_n.*L./v_n)*F_n*c_n./(F_n*c_n+K_n)),[0 1]);
        F_nall = [F_nall;F_n];
        c_r = rho*Y(k,ntot+2)/((1-1/a)*M);
        F_r = fzero(@(F_r) Y(k,ntot+2)*(1-F_r)-sum(Y(k,1:ntot)'.*(1+Gamma_r.*L./v_r)*F_r*c_r./(F_r*c_r+K_r)),[0 1]);
        F_rall = [F_rall;F_r];
    end
    Y0 = Y(end,:)'*gamma;
    Y0(ntot+4) = Y(end,ntot+4)*0.45;
    tlast = t_all(end);
    g = g/2;
    tend(generation) = length(t_all); % record the timepoint at division
end