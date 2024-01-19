% This user defined function will give the surface sag, Z(h) where
% h=sqrt(x^2 +y^2), given KrElt, KcElt, x-coord, and y-coord. The sag
% equation for conic surfaces defines all  points on the surface in 3-D
% space
%
% Variables:
%
% KrElt   -radius of curvature of conic
% KcElt   -conic constant KcElt<-1 hyperboloid, KcElt=-1 paraboloid
%           -1 < KcElt <0 ellipsoid with major axis on the principal axis
%           KcElt=0 sphere, and KcElt>0 oblate ellipsoid 
%         
% x       -x-axis coordinate
% y       -y-axis coordinate
% h       -radial pupil coordinate, the magnitude projection of the point
%          of incidence vector onto the plane perpendicular to psiElt
function[Z]=AsphericSag(KrElt,KcElt,a1,a2,a3,x,y)
r=sqrt(x.^2 + y.^2);
Z=((r.^2)./(KrElt+sqrt(KrElt.^2 -((1+KcElt).*r^.2))))+a1*r.^4+a2*r.^6+a3*r.^8;
return