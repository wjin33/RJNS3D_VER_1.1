function [r_v_m,hcuboid]=GenCuboid(x0,y0,z0,dx,dy,dz)
%(x0,y0,z0)--coordinate of the first vertex
%(dx,dy,dz)--length, width, and height og the cuboid
%r_v_m--24x5 Matrix, 24x(1~3) are the coordinates of each face, 24x(4~5) are the attitude of each face
hcuboid=figure;    %handle of the cuboid figure
colort=[0,0,1];  %color of the cuboid figure

%subsurface£¬generate vertexes from the first, anticlockwise
x=[x0,x0+dx,x0+dx,x0+dx,x0+dx,x0,   x0,   x0];
y=[y0,y0,   y0,   y0+dy,y0+dy,y0+dy,y0+dy,y0];
z=[z0,z0,   z0,   z0,   z0,   z0,   z0,   z0];
line(x,y,z,'Color',colort)
hold on

%topsurface£¬generate vertexes from the first, anticlockwise
x=[x0,x0+dx,x0+dx,x0+dx,x0+dx,x0,x0,x0];
y=[y0,y0,y0,y0+dy,y0+dy,y0+dy,y0+dy,y0];
z=[z0+dz,z0+dz,z0+dz,z0+dz,z0+dz,z0+dz,z0+dz,z0+dz];
line(x,y,z,'Color',colort)
hold on

%generate ridges from the first point, anticlockwise
x=[x0,x0];
y=[y0,y0];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

x=[x0+dx,x0+dx];
y=[y0,y0];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

x=[x0+dx,x0+dx];
y=[y0+dy,y0+dy];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

x=[x0,x0];
y=[y0+dy,y0+dy];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

%initialize the Matrix which stored the coordinates and attitude of each surface in row 
r_v_m=zeros(24,5);

%face 1, the coordinates and attitude of topsurface 
x=x0; y=y0; z=z0+dz;
r_v_m(1,1:3)=[x,y,z];
x=x0+dx; y=y0; z=z0+dz;
r_v_m(2,1:3)=[x,y,z];
x=x0+dx; y=y0+dy; z=z0+dz;
r_v_m(3,1:3)=[x,y,z];
x=x0; y=y0+dy; z=z0+dz;
r_v_m(4,1:3)=[x,y,z];
%dip and trend
for i=1:4
	r_v_m(i,4:5)=[pi/2.0,0];
end

%face 2
x=x0; y=y0+dy; z=z0+dz;
r_v_m(5,1:3)=[x,y,z];
x=x0; y=y0+dy; z=z0;
r_v_m(6,1:3)=[x,y,z];
x=x0; y=y0; z=z0;
r_v_m(7,1:3)=[x,y,z];
x=x0; y=y0; z=z0+dz;
r_v_m(8,1:3)=[x,y,z];
for i=5:8
	r_v_m(i,4:5)=[1.5*pi,pi/2.0];
end

%face 3
x=x0; y=y0; z=z0+dz;
r_v_m(9,1:3)=[x,y,z];
x=x0; y=y0; z=z0;
r_v_m(10,1:3)=[x,y,z];
x=x0+dx; y=y0; z=z0;
r_v_m(11,1:3)=[x,y,z];
x=x0+dx; y=y0; z=z0+dz;
r_v_m(12,1:3)=[x,y,z];
for i=9:12
	r_v_m(i,4:5)=[pi,pi/2.0];
end

%face 4
x=x0+dx; y=y0; z=z0+dz;
r_v_m(13,1:3)=[x,y,z];
x=x0+dx; y=y0; z=z0;
r_v_m(14,1:3)=[x,y,z];
x=x0+dx; y=y0+dy; z=z0;
r_v_m(15,1:3)=[x,y,z];
x=x0+dx; y=y0+dy; z=z0+dz;
r_v_m(16,1:3)=[x,y,z];
for i=13:16
	r_v_m(i,4:5)=[pi/2.0,pi/2.0];
end

%face 5
x=x0+dx; y=y0+dy; z=z0+dz;
r_v_m(17,1:3)=[x,y,z];
x=x0+dx; y=y0+dy; z=z0;
r_v_m(18,1:3)=[x,y,z];
x=x0; y=y0+dy; z=z0;
r_v_m(19,1:3)=[x,y,z];
x=x0; y=y0+dy; z=z0+dz;
r_v_m(20,1:3)=[x,y,z];
for i=17:20
	r_v_m(i,4:5)=[0,pi/2.0];
end

%face 6, subsurface
x=x0+dx; y=y0; z=z0;
r_v_m(21,1:3)=[x,y,z];
x=x0; y=y0; z=z0;
r_v_m(22,1:3)=[x,y,z];
x=x0; y=y0+dy; z=z0;
r_v_m(23,1:3)=[x,y,z];
x=x0+dx; y=y0+dy; z=z0;
r_v_m(24,1:3)=[x,y,z];
for i=21:24
	r_v_m(i,4:5)=[3.0*pi/2.0,0];
end