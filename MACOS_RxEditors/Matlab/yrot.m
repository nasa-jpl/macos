%Given angle betta it generates the y-axis rotation matrix yrotm
function[yrotm]=yrot(betta)
yrotm=[cos(betta) 0 sin(betta);0 1 0; -sin(betta) 0 cos(betta)];
return