function GenCuboidDis(sx0,sy0,sz0,si,sj,sk,sx,sy,sz,colort)
%(sx,sy,sz)--length,width,hight of the cuboid
%sx0,sy0,sz0--coordinate of cuboid
%si,sj,sk--initial offset

x0=sx0+(si-1)*sx;  y0=sy0+(sj-1)*sy;  z0=sz0+(sk-1)*sz;
%subsurface
x=[x0, x0+sx, x0+sx, x0+sx, x0+sx, x0,    x0,    x0];
y=[y0, y0,    y0,    y0+sy, y0+sy, y0+sy, y0+sy, y0];
z=[z0, z0,    z0,    z0,    z0,    z0,    z0,    z0];
line(x,y,z,'Color',colort,'LineWidth',2)
hold on

%top surface
x=[x0,   x0+sx, x0+sx, x0+sx, x0+sx, x0,    x0,      x0];
y=[y0,   y0,    y0,    y0+sy, y0+sy, y0+sy, y0+sy,   y0];
z=[z0+sz,z0+sz, z0+sz, z0+sz, z0+sz, z0+sz, z0+sz,   z0+sz];
line(x,y,z,'Color',colort,'LineWidth',2)
hold on

%generate the ridges, anticlockwise
x=[x0,x0];
y=[y0,y0];
z=[z0,z0+sz];
line(x,y,z,'Color',colort,'LineWidth',2)
hold on

x=[x0+sx,x0+sx];
y=[y0,y0];
z=[z0,z0+sz];
line(x,y,z,'Color',colort,'LineWidth',2)
hold on

x=[x0+sx,x0+sx];
y=[y0+sy,y0+sy];
z=[z0,z0+sz];
line(x,y,z,'Color',colort,'LineWidth',2)
hold on

x=[x0,x0];
y=[y0+sy,y0+sy];
z=[z0,z0+sz];
line(x,y,z,'Color',colort,'LineWidth',2)
