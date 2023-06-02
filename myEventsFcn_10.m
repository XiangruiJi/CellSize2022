function [position,isterminal,direction] = myEventsFcn_10(t,Y,ntot,rho,L,threshold)
% inhibitor-dilution
position = (threshold-rho*Y(ntot+4)/sum(Y(ntot+1:2*ntot).*L))+(t<5); % The value that we want to be zero.
isterminal = 1;  % Halt integration
direction = 1;   % The zero can only be approached from one direction
end