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
% h       -radial distance to projection  of the point of incidence onto a
%          plane perpendicular to the principal axis and including surface
%          vertex
function[Z]=ConicSag(KrElt,KcElt,x,y)
h=sqrt(x.^2 + y.^2);
Z=((h.^2)./(KrElt+sqrt(KrElt.^2 -((1+KcElt).*(h.^2)))));
return