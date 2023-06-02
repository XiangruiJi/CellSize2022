function [position,isterminal,direction] = myEventsFcn_12(t,Y,ntot,rho,L,thresholdact,thresholdinh)
% AND gate, one activator and one inhibitor
position = min([rho*Y(ntot+3)/sum(Y(ntot+1:2*ntot).*L)-thresholdact;thresholdinh-rho*Y(ntot+4)/sum(Y(ntot+1:2*ntot).*L)])+(t<5); % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 1;   % The zero can only be approached from one direction
end