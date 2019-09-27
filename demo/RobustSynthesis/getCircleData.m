function [X,Y] = getCircleData(c,r)
% Creates Circle Data X and Y coordinates from given center and radius
theta = linspace(-pi,pi,1000);
c1 = c(1);
c2 = c(2);
X = c1 + r*cos(theta);
Y = c2 + r*sin(theta);
end