%Given angle alpha in Radians it generates the x-axis rotation matrix xrotm
function[xrotm]=xrot(alpha)
xrotm=[1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)]; %rotation about x-axis in the clockwise direction
%xrotm=[1 0 0; cos(alpha) sin(alpha) 0;-sin(alpha) cos(alpha) 0]; %changed on April 30th, 2010
return
