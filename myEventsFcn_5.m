function [position,isterminal,direction] = myEventsFcn_5(t,Y,ntot,threshold,gx)
position = (sum(Y(ntot+3:ntot+2+gx))/sum(Y(ntot+3+gx:ntot+2+2*gx))-threshold)+(t<5); % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 1;   % The zero can only be approached from one direction
end