%Given angle gamma in Radians it generates the z-axis rotation matrix zrotm

function[zrotm]=zrot(gamma)
zrotm=[cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
return