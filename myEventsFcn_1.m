function [position,isterminal,direction] = myEventsFcn_1(t,Y,Mmax,ntot,L)
position = sum(Y(ntot+1:2*ntot).*L)-Mmax; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % The zero can only be approached from one direction
end