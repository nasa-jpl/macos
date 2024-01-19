function[Z]=AsphSurf(KrElt,KcElt,a1,a2,a3,x,y)
r=sqrt(x.^2 + y.^2);
Z=((r.^2)./(KrElt+sqrt(KrElt.^2 -((1+KcElt).*r^.2))))+a1*r.^4+a2*r.^6+a3*r.^8;
return