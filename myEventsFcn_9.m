function [position,isterminal,direction] = myEventsFcn_9(t,Y,ntot,rho,L,threshold)
% activator-acummulation
position = (rho*Y(ntot+3)/sum(Y(ntot+1:2*ntot).*L)-threshold)+(t<5); % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 1;   % The zero can only be approached from one direction
end