% Octavio Narváez Aroche
% Berkeley Center for Control and Identification
% Summer 2017
% Events function for the TV LQR control of a three link robot modeling the 
% STS movement of an exoskeleton and its user.

function [event,isterminal,direction] = TVLQRSTS3LinkEvents(t,x,xmin,xmax)

% Number of states 
nx = numel(xmin);

% "Zero" events.
event(1:nx) = x(1:nx)-xmin;
event(nx+1:2*nx) = xmax-x(1:nx); 

% Halt integration when events are triggered.
isterminal(1:2*nx) = 1;

% The zero can be approached from above (event is decreasing to 0).
direction = -1;