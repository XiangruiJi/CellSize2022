function [position,isterminal,direction] = myEventsFcn_7(t,Y,ntot,threshold)
position = Y(ntot+3)/Y(ntot+4)-threshold; % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 1;   % The zero can only be approached from one direction
end