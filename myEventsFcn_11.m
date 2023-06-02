function [position,isterminal,direction] = myEventsFcn_11(t,Y,ntot,threshold,gx)
position = ((2*Y(ntot+3)+2*Y(ntot+4)+sum(Y(ntot+5:ntot+gx-2))+0.5*sum(Y(ntot+gx-1:ntot+2+gx)))/(2*Y(ntot+3+gx)+2*Y(ntot+4+gx)+sum(Y(ntot+5+gx:ntot+2*gx-2))+0.5*sum(Y(ntot+2*gx-1:ntot+2+2*gx)))-threshold)+(t<5); % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 1;   % The zero can only be approached from one direction
end