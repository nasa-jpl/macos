% Applies a rotation theta in Radians of the vector X about the axis or 
% vector given by n.
% This performs an axis angle rotation given the angle theta of rotation,
% the vector containing the point x to be rotated, and n which is the axis 
% of rotation of the point x

function[newpoint]=axisanglerot(x,theta,n)
newpoint=x.*cos(theta)+n*(dot(n,x))*(1-cos(theta))+(cross(x,n)*sin(theta));
return





