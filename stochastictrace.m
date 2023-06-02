function [meanV,CVV] = stochastictrace(thetabar,dd,genmax,m_0,p_0,tmax,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p,gx)
meanV = NaN;
CVV = NaN;
% simulate the cell cycle to a steady state using the deterministic model
options = odeset('Events',@(t,Y) myEventsFcn_5(t,Y,ntot,thetabar,gx));
Y0 = [m_0;p_0];
for generation = 1:15
    [t,Y,~,~,~] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,tmax],Y0,options);
    Y0 = Y(end,:)'/2;
    if t(end)>=tmax
        break;
    end
end
if generation==15
    % the stochastic model
    Vb = NaN(genmax,1);
    theta = thetabar+dd/3*randn(genmax,1); % \theta~N(theta,(dd/3)^2)
    for generation = 1:genmax
        if theta(generation)<thetabar-dd || theta(generation)>thetabar+dd
            theta(generation) = thetabar;
        end
    end
    for generation = 1:genmax
        options = odeset('Events',@(t,Y) myEventsFcn_5(t,Y,ntot,theta(generation),gx));
        [t,Y,~,~,~] = ode45(@(t,Y) odefun(t,Y,ntot,rho,a,g,L,Gamma_n,K_n,v_n,tau_m,Gamma_r,K_r,v_r,tau_p),[0,tmax],Y0,options);
        Y0 = Y(end,:)'/2;
        Vb(generation) = sum(Y(1,ntot+1:2*ntot)'.*L)/rho;
        if t(end)>=tmax
            break;
        end
    end
    if generation==genmax
        meanV = mean(Vb);
        CVV = std(Vb)/mean(Vb);
    end
end
end