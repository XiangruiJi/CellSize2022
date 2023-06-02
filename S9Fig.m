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
K_n = 6000*ones(ntot,1); % (um^(-3))
K_n(3) = 12000;
K_n(4) = 4000;
v_n = 720; % (aa/min)
tau_m = 5*ones(ntot,1); % (min)
Gamma_r = 10*ones(ntot,1); % (min^(-1))
K_r = 6000*ones(ntot,1); % (um^(-3))
v_r = 720; % (aa/min)
tau_p = 1e4*ones(ntot,1); % (min)
tau_p(3:4) = 5;

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

genmax = 5000; % the number of generations
theta1 = K_n(4)/K_n(3);
theta2 = (c_n+K_n(4))/(c_n+K_n(3));
C1 = 200;
C2 = 1;
% First, simulate the cell cycle to a steady state using the deterministic model
options = odeset('Events',@(t,Y) myEventsFcn_3(t,Y,ntot,0.7));
Y0 = [m_0;p_0];
for generation = 1:10
    [t,Y,~,~,~] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,tmax],Y0,options);
    if t(end)>=tmax
        break;
    end
    Y0 = Y(end,:)'/2;
end
if generation==10
    % the stochastic model
    V_birth = zeros(genmax,1);
    Vd = zeros(genmax,1);
    thetapre = rand(genmax,1); % \thetapre~U(0,1)
    for generation = 1:genmax
        [~,Y] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,20],Y0);
        V_birth(generation) = sum(Y(1,ntot+1:2*ntot)'.*L)/rho;
        Y0 = Y(end,:)';
        thetab = Y0(ntot+3)/Y0(ntot+4);
        theta = fzero(@(x) log(1-thetapre(generation))+integral(@(x1) C1*((x1-theta1)./(theta2-x1)).^C2,thetab,x),[thetab,theta2]); % inverse transform sampling
        options = odeset('Events',@(t,Y) myEventsFcn_7(t,Y,ntot,theta));
        [t,Y,~,~,~] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,tmax],Y0,options);
        if t(end)>=tmax
            break;
        end
        Vd(generation) = sum(Y(end,ntot+1:2*ntot)'.*L)/rho;
        Y0 = Y(end,:)'/2;
    end
end