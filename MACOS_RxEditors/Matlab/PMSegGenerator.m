%==========================================================================
%==========================================================================
% Purpose: The script purpose is to calculate an approximated segmented 
%          mirror in place of a monolithic primary mirror. We make use of
%          the monolithic primary parameters Kr, Kc, and VptElt, along with
%          the segment size information to perform the approximation. The
%          script generates segment vertices, segment center, and segment
%          TElt better known as the global to local transformation matrix
%
% Record:
%     Date                  Author               Description
% ===============      ==================    ==============================
%  1/23/2013           Luis F. Marchen       Modified Version of LD's
%                                            Segment Generator
%
% Define Variables and Constants:
%
% Hd           -- vertext to vertex distance, units required meters
% He           -- side to opposing side distance, units required meters
% gap          -- size of gab between segments, units required meters
% Kr           -- radius of curvature, units required meters
% Kc           -- conic constant
% ifSurface    -- specifies the type of primary mirror surface Conic or 
%                 Aspheric
%
%==========================================================================
%==========================================================================
function[seg_data,SegCen,TElt]=PMSegGenerator(Hd,He,Kc,Kr,gap,ifSurface)
%--------------------------------------------------------------------------
% Primary mirror surface and segment parameters
%--------------------------------------------------------------------------
h=Hd/2;       % hexagon radius, diamter=point-to-point distance of hexagon (m)
R=Kr;   % radius of curvature (m)
hd=Hd/2;
he=He/2;
rhoA=gap+2*he;
rhoB=2*gap+4*he;
rhoC=((3*gap/2)+3*he)/(cos(pi/6));
%rho=[1.234134,2.4682686,2.137583];%segment offset distances for A, B, and C segs, contains gap size with in it
rho=[rhoA, rhoB, rhoC];
group=['A', 'B', 'C'];

%--------------------------------------------------------------------------
% define arrays to save vertex data for segments
%--------------------------------------------------------------------------
SegAVtxMat=zeros(36,6);
SegBVtxMat=zeros(36,6);
SegCVtxMat=zeros(36,6);
%--------------------------------------------------------------------------
% define aspheric coefficients here if the PM surface shape is an aspheric
% surface. The user defined function to compute aspheric surface sag is
% only set up to handle three aspheric coefficients, if more are needed
% user must edit AsphericSag.m
%--------------------------------------------------------------------------
%a1=;
%a2=;
%a3=;
%==========================================================================
% open output file
%==========================================================================
timdatstr = date;
filename=strcat('segment',timdatstr);
fid=fopen(filename,'wt');
fprintf(fid,'%s \n',filename);
fprintf(fid,'hexagon point-to-point diameter (m) %11.7f\n',h);
fprintf(fid,'parabola radius of curvature (m) %11.7f\n',R);
fprintf(fid,'segment offsets (m) %11.7f %11.7f %11.7f \n',rho);
figure(1)
figure(2)
%==========================================================================
%==========================================================================
% generate segment vertices and centers
%==========================================================================
%==========================================================================
for i=1:3
%==========================================================================
% A and B segment vertices
%==========================================================================
    if (i==1 || i==2)
        fprintf(fid,'%c segments. Offset = %11.7f (m) \n', group(i),rho(i));
        fprintf(fid,'phi x y z \n');
        %for phi=0.:pi/3.:5.*pi/3.; %[2*pi/3, pi/3, 0, 5*pi/3, 4*pi/3 3*pi/3]
        vtxite=0; %veterx iteration for each segment
        for phi=[2*pi/3, pi/3, 0, 5*pi/3, 4*pi/3 3*pi/3];
            vtxite=vtxite+1;
%--------------------------------------------------------------------------
%   coordinates of points for standard hexagon
%--------------------------------------------------------------------------
            deg30=pi/6;
            h0=[0, 0, 0]';
            h1=[h, 0, 0]';
            h2=[h*sin(deg30), h*cos(deg30), 0]';
            h3=[-h*sin(deg30), h*cos(deg30), 0]';
            h4=[-h, 0, 0]';
            h5=[-h*sin(deg30), -h*cos(deg30), 0]';
            h6=[h*sin(deg30), -h*cos(deg30), 0]';          
%--------------------------------------------------------------------------
% rotate about panel x-axis by atan(rho/R)
%--------------------------------------------------------------------------
            psi=atan(rho(i)/R);
            R1=[1, 0, 0;0 cos(psi), -sin(psi);0, sin(psi), cos(psi)];
            hh0=R1*h0;
            hh1=R1*h1;
            hh2=R1*h2;
            hh3=R1*h3;
            hh4=R1*h4;
            hh5=R1*h5;
            hh6=R1*h6;          
%--------------------------------------------------------------------------
% translate center to (-rho,phi=0)
%--------------------------------------------------------------------------
            delta=[0., -rho(i), 0]';
            hhh0=hh0+delta;
            hhh1=hh1+delta;
            hhh2=hh2+delta;
            hhh3=hh3+delta;
            hhh4=hh4+delta;
            hhh5=hh5+delta;
            hhh6=hh6+delta;
%--------------------------------------------------------------------------
% rotate about z-axis by angle phi
%--------------------------------------------------------------------------
            R2=[cos(phi), sin(phi), 0;-sin(phi), cos(phi), 0;0, 0, 1];
            hhhh0=R2*hhh0;
            hhhh1=R2*hhh1;
            hhhh2=R2*hhh2;
            hhhh3=R2*hhh3;
            hhhh4=R2*hhh4;
            hhhh5=R2*hhh5;
            hhhh6=R2*hhh6;
%--------------------------------------------------------------------------
% Compute local axes for A and B segments
%--------------------------------------------------------------------------
            ax=[1 0 0];
            ay=[0 1 0];
            az=[0 0 1];
            
            aax=R1*ax';
            aay=R1*ay';
            aaz=R1*az';
            
            aaax=R2*aax;
            aaay=R2*aay;
            aaaz=R2*aaz;           
%--------------------------------------------------------------------------
% The surface sag is computed below based on surface type: it can be a
% conic surface or an aspheric surface
%--------------------------------------------------------------------------
        if strcmp(ifSurface,'Conic');
            z0=-ConicSag(Kr,Kc,hhhh0(1),hhhh0(2));
            z1=-ConicSag(Kr,Kc,hhhh1(1),hhhh1(2));
            z2=-ConicSag(Kr,Kc,hhhh2(1),hhhh2(2));
            z3=-ConicSag(Kr,Kc,hhhh3(1),hhhh3(2));
            z4=-ConicSag(Kr,Kc,hhhh4(1),hhhh4(2));
            z5=-ConicSag(Kr,Kc,hhhh5(1),hhhh5(2));
            z6=-ConicSag(Kr,Kc,hhhh6(1),hhhh6(2));
            
        elseif strcmpi(ifSurface,'Aspheric')
            z0=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh0(1),hhhh0(2));
            z1=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh1(1),hhhh1(2));
            z2=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh2(1),hhhh2(2));
            z3=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh3(1),hhhh3(2));
            z4=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh4(1),hhhh4(2));
            z5=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh5(1),hhhh5(2));
            z6=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh6(1),hhhh6(2));
        end
%--------------------------------------------------------------------------
% form arrays for plotting
%--------------------------------------------------------------------------
            x=[hhhh1(1), hhhh2(1), hhhh3(1), hhhh4(1), hhhh5(1), hhhh6(1), hhhh1(1)];
            y=[hhhh1(2), hhhh2(2), hhhh3(2), hhhh4(2), hhhh5(2), hhhh6(2), hhhh1(2)];
            z=[z1, z2, z3, z4, z5, z6, z1];
            zz=[hhhh1(3), hhhh2(3), hhhh3(3), hhhh4(3), hhhh5(3), hhhh6(3), hhhh1(3)];
            if (i==1)
                figure(1)
                plot3(x,y,z,'r')
                hold on
                plot3(hhhh0(1),hhhh0(2),z0,'r*')
                figure(2)
                plot3(x,y,zz,'r')
                hold on
            else
                figure(1)
                plot3(x,y,z,'g')
                plot3(hhhh0(1),hhhh0(2),z0,'g*')
                figure(2)
                plot3(x,y,zz,'g')
            end
%--------------------------------------------------------------------------
% compose matrix of vertecies for each segement A and B
%--------------------------------------------------------------------------           
            if i==1
                SegAVtxMat((vtxite-1)*6+1:6*vtxite,1:6)=[hhhh3(1),hhhh3(2),hhhh3(3), hhhh3(1),hhhh3(2),z3;...
                                                         hhhh4(1),hhhh4(2),hhhh4(3), hhhh4(1),hhhh4(2),z4;...
                                                         hhhh5(1),hhhh5(2),hhhh5(3), hhhh5(1),hhhh5(2),z5;...
                                                         hhhh6(1),hhhh6(2),hhhh6(3), hhhh6(1),hhhh6(2),z6;...
                                                         hhhh1(1),hhhh1(2),hhhh1(3), hhhh1(1),hhhh1(2),z1;...
                                                         hhhh2(1),hhhh2(2),hhhh2(3), hhhh2(1),hhhh2(2),z2];
                
                SegACenMat(vtxite,:)=[hhhh0(1),hhhh0(2),z0];   
                TEltA(:,:,vtxite)=[[aaax;zeros(3,1)],[aaay;zeros(3,1)],[aaaz;zeros(3,1)],[zeros(3,1);aaax],[zeros(3,1);aaay],[zeros(3,1);aaaz]];
                                                     
            elseif i==2
                 SegBVtxMat((vtxite-1)*6+1:6*vtxite,1:6)=[hhhh3(1),hhhh3(2),hhhh3(3), hhhh3(1),hhhh3(2),z3;...
                                                          hhhh4(1),hhhh4(2),hhhh4(3), hhhh4(1),hhhh4(2),z4;...
                                                          hhhh5(1),hhhh5(2),hhhh5(3), hhhh5(1),hhhh5(2),z5;...
                                                          hhhh6(1),hhhh6(2),hhhh6(3), hhhh6(1),hhhh6(2),z6;...
                                                          hhhh1(1),hhhh1(2),hhhh1(3), hhhh1(1),hhhh1(2),z1;...
                                                          hhhh2(1),hhhh2(2),hhhh2(3), hhhh2(1),hhhh2(2),z2];
                                                      
                 SegBCenMat(vtxite,:)=[hhhh0(1),hhhh0(2),z0];
                 TEltB(:,:,vtxite)=[[aaax;zeros(3,1)],[aaay;zeros(3,1)],[aaaz;zeros(3,1)],[zeros(3,1);aaax],[zeros(3,1);aaay],[zeros(3,1);aaaz]];
            end
%--------------------------------------------------------------------------
% print results
%--------------------------------------------------------------------------
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh3(1),hhhh3(2),z3);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh4(1),hhhh4(2),z4);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh5(1),hhhh5(2),z5);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh6(1),hhhh6(2),z6);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh1(1),hhhh1(2),z1);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh2(1),hhhh2(2),z2);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f center \n',phi*180/pi,hhhh0(1),hhhh0(2),z0);
        end
    else
%==========================================================================
% generate C segment vertices
%==========================================================================
        fprintf(fid,'%c segments. Offset = %11.7f (m) \n', group(3),rho(3));
        fprintf(fid,'phi x y z \n');
        vtxite=0;
        %for phi=pi/6:pi/3.:11.*pi/6.;
        for phi=[pi/6 + pi/3, pi/6, pi/6 + 5*pi/3, pi/6 + 4*pi/3, pi/6 + pi,pi/6 + 2*pi/3 ];
            
            vtxite=vtxite+1;
%--------------------------------------------------------------------------
%   coordinates of points for standard hexagon
%--------------------------------------------------------------------------
            deg30=pi/6;
            h0=[0, 0, 0]';
            h1=[h, 0, 0]';
            h2=[h*sin(deg30), h*cos(deg30), 0]';
            h3=[-h*sin(deg30), h*cos(deg30), 0]';
            h4=[-h, 0, 0]';
            h5=[-h*sin(deg30), -h*cos(deg30), 0]';
            h6=[h*sin(deg30), -h*cos(deg30), 0]';
%--------------------------------------------------------------------------
% rotate about z-axis by 30 degrees to get body x-axis pointing at a
% flat and not a point
%--------------------------------------------------------------------------
            R0=[cos(deg30), sin(deg30), 0;-sin(deg30), cos(deg30), 0;0, 0, 1];
            h0=R0*h0;
            h1=R0*h1;
            h2=R0*h2;
            h3=R0*h3;
            h4=R0*h4;
            h5=R0*h5;
            h6=R0*h6;
%--------------------------------------------------------------------------
% rotate about panel x-axis by atan(rho/R)
%--------------------------------------------------------------------------
            psi=atan(rho(3)/R);
            R1=[1, 0, 0;0 cos(psi), -sin(psi);0, sin(psi), cos(psi)];
            hh0=R1*h0;
            hh1=R1*h1;
            hh2=R1*h2;
            hh3=R1*h3;
            hh4=R1*h4;
            hh5=R1*h5;
            hh6=R1*h6;
%--------------------------------------------------------------------------
% translate center to (rho,phi=0)
%--------------------------------------------------------------------------
            delta=[0, -rho(3), 0]';
            hhh0=hh0+delta;
            hhh1=hh1+delta;
            hhh2=hh2+delta;
            hhh3=hh3+delta;
            hhh4=hh4+delta;
            hhh5=hh5+delta;
            hhh6=hh6+delta;
%--------------------------------------------------------------------------
% rotate about z-axis by angle phi
%--------------------------------------------------------------------------
            R2=[cos(phi), sin(phi), 0;-sin(phi), cos(phi), 0;0, 0, 1];
            hhhh0=R2*hhh0;
            hhhh1=R2*hhh1;
            hhhh2=R2*hhh2;
            hhhh3=R2*hhh3;
            hhhh4=R2*hhh4;
            hhhh5=R2*hhh5;
            hhhh6=R2*hhh6;
%--------------------------------------------------------------------------
% Compute local axes for A and B segments
%--------------------------------------------------------------------------
            ax=[1 0 0];
            ay=[0 1 0];
            az=[0 0 1];
            
            aax=R1*ax';
            aay=R1*ay';
            aaz=R1*az';
            
            aaax=R2*aax;
            aaay=R2*aay;
            aaaz=R2*aaz;
            
            TEltC(:,:,vtxite)=[[aaax;zeros(3,1)],[aaay;zeros(3,1)],[aaaz;zeros(3,1)],[zeros(3,1);aaax],[zeros(3,1);aaay],[zeros(3,1);aaaz]];
%--------------------------------------------------------------------------
% calculate the height of the paraboloid at each of these points
%--------------------------------------------------------------------------
         if strcmp(ifSurface,'Conic');
            z0=-ConicSag(Kr,Kc,hhhh0(1),hhhh0(2));
            z1=-ConicSag(Kr,Kc,hhhh1(1),hhhh1(2));
            z2=-ConicSag(Kr,Kc,hhhh2(1),hhhh2(2));
            z3=-ConicSag(Kr,Kc,hhhh3(1),hhhh3(2));
            z4=-ConicSag(Kr,Kc,hhhh4(1),hhhh4(2));
            z5=-ConicSag(Kr,Kc,hhhh5(1),hhhh5(2));
            z6=-ConicSag(Kr,Kc,hhhh6(1),hhhh6(2));
            
        elseif strcmpi(ifSurface,'Aspheric')
            z0=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh0(1),hhhh0(2));
            z1=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh1(1),hhhh1(2));
            z2=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh2(1),hhhh2(2));
            z3=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh3(1),hhhh3(2));
            z4=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh4(1),hhhh4(2));
            z5=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh5(1),hhhh5(2));
            z6=-AsphSurf(Kr,Kc,a1,a2,a3,hhhh6(1),hhhh6(2));
         end
%--------------------------------------------------------------------------
% form arrays for plotting
%--------------------------------------------------------------------------
            x=[hhhh1(1), hhhh2(1), hhhh3(1), hhhh4(1), hhhh5(1), hhhh6(1), hhhh1(1)];
            y=[hhhh1(2), hhhh2(2), hhhh3(2), hhhh4(2), hhhh5(2), hhhh6(2), hhhh1(2)];
            z=[z1, z2, z3, z4, z5, z6, z1];
            zz=[hhhh1(3), hhhh2(3), hhhh3(3), hhhh4(3), hhhh5(3), hhhh6(3), hhhh1(3)];
            figure(1)
            plot3(x,y,z)
            plot3(hhhh0(1),hhhh0(2),z0,'b*')
            figure(2)
            plot3(x,y,zz)          
%--------------------------------------------------------------------------
% compose matrix of vertecies for C segements
%--------------------------------------------------------------------------
            SegCVtxMat((vtxite-1)*6+1:6*vtxite,1:6)=[hhhh3(1), hhhh3(2), hhhh3(3), hhhh3(1),hhhh3(2),z3;...
                                                     hhhh4(1), hhhh4(2), hhhh4(3), hhhh4(1),hhhh4(2),z4;...
                                                     hhhh5(1), hhhh5(2), hhhh5(3), hhhh5(1),hhhh5(2),z5;...
                                                     hhhh6(1), hhhh6(2), hhhh6(3), hhhh6(1),hhhh6(2),z6;...
                                                     hhhh1(1), hhhh1(2), hhhh1(3), hhhh1(1),hhhh1(2),z1;...
                                                     hhhh2(1), hhhh2(2), hhhh2(3), hhhh2(1),hhhh2(2),z2];
                                                         
           SegCCenMat(vtxite,:)=[hhhh0(1),hhhh0(2),z0];                                                         
%--------------------------------------------------------------------------
% print results
%--------------------------------------------------------------------------
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh3(1),hhhh3(2),z3);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh4(1),hhhh4(2),z4);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh5(1),hhhh5(2),z5);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh6(1),hhhh6(2),z6);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh1(1),hhhh1(2),z1);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f \n',phi*180/pi,hhhh2(1),hhhh2(2),z2);
            fprintf(fid,'%5.1f %11.14f %11.14f %11.14f center \n',phi*180/pi,hhhh0(1),hhhh0(2),z0);

        end
    end
end
%==========================================================================
% compose Segment_Vertices matrix
%==========================================================================
SegBCVtx=zeros(72,6);
for jj=1:6;
    SegBCVtx((1:6)+(12*(jj-1)),1:6)=SegBVtxMat((jj-1)*6+1:6*jj,1:6);
    SegBCVtx((7:12)+(12*(jj-1)),1:6)=SegCVtxMat((jj-1)*6+1:6*jj,1:6);
end

SegVtx=[SegAVtxMat;SegBCVtx]; %segment vertices
seg_vertices=SegVtx(:,4:6);
%==========================================================================
% TElt matrix
%==========================================================================
TElt=zeros(6,6,18);
SegCen=zeros(18,3);

tite=0;
for uu=1:6;
    tite=tite+1;
    TElt(:,:,uu)=TEltA(:,:,tite); 
    SegCen(uu ,:)=SegACenMat(tite,:);
end

tite=0;
for uu=[7 9 11 13 15 17];
    tite=tite+1;
    TElt(:,:,uu)=TEltB(:,:,tite); 
    SegCen(uu ,:)=SegBCenMat(tite,:);
end

tite=0;
for uu=[8 10 12 14 16 18];
    tite=tite+1;
    TElt(:,:,uu)=TEltC(:,:,tite); 
    SegCen(uu ,:)=SegCCenMat(tite,:);
end

seg_data=zeros(6,3,18);
for hh=1:18;
    seg_data(:,:,hh)=seg_vertices( (1+6*(hh-1)):6*(hh),:);
end
%==========================================================================
% segment vertex labels
%==========================================================================
figure(1)
grid
xlabel('x')
ylabel('y')
zlabel('z')
figure(2)
grid
xlabel('x')
ylabel('y')
zlabel('z')
fclose(fid);

SegVtxLabel={'A1_1' ;   'A1_2' ;   'A1_3'  ;  'A1_4'  ;  'A1_5'  ;  'A1_6';...
     'B1_1' ;   'B1_2' ;   'B1_3'  ;  'B1_4'  ;  'B1_5'  ;  'B1_6';...
     'C1_1' ;   'C1_2' ;   'C1_3'  ;  'C1_4'  ;  'C1_5'  ;  'C1_6';...
     'A2_1' ;   'A2_2' ;   'A2_3'  ;  'A2_4'  ;  'A2_5'  ;  'A2_6';...
     'B2_1' ;   'B2_2' ;   'B2_3'  ;  'B2_4'  ;  'B2_5'  ;  'B2_6';...
     'C2_1' ;   'C2_2' ;   'C2_3'  ;  'C2_4'  ;  'C2_5'  ;  'C2_6';...
     'A3_1' ;   'A3_2' ;   'A3_3'  ;  'A3_4'  ;  'A3_5'  ;  'A3_6';...
     'B3_1' ;   'B3_2' ;   'B3_3'  ;  'B3_4'  ;  'B3_5'  ;  'B3_6';...
     'C3_1' ;   'C3_2' ;   'C3_3'  ;  'C3_4'  ;  'C3_5'  ;  'C3_6';...
     'A4_1' ;   'A4_2' ;   'A4_3'  ;  'A4_4'  ;  'A4_5'  ;  'A4_6';...
     'B4_1' ;   'B4_2' ;   'B4_3'  ;  'B4_4'  ;  'B4_5'  ;  'B4_6';...
     'C4_1' ;   'C4_2' ;   'C4_3'  ;  'C4_4'  ;  'C4_5'  ;  'C4_6';...
     'A5_1' ;   'A5_2' ;   'A5_3'  ;  'A5_4'  ;  'A5_5'  ;  'A5_6';...
     'B5_1' ;   'B5_2' ;   'B5_3'  ;  'B5_4'  ;  'B5_5'  ;  'B5_6';...
     'C5_1' ;   'C5_2' ;   'C5_3'  ;  'C5_4'  ;  'C5_5'  ;  'C5_6';...
     'A6_1' ;   'A6_2' ;   'A6_3'  ;  'A6_4'  ;  'A6_5'  ;  'A6_6';...
     'B6_1' ;   'B6_2' ;   'B6_3'  ;  'B6_4'  ;  'B6_5'  ;  'B6_6';...
     'C6_1' ;   'C6_2' ;   'C6_3'  ;  'C6_4'  ;  'C6_5'  ;  'C6_6'};
return
