function Genborehole(location,bore_dip,H,R,m)
%绘制任意位置、大小,产状的圆柱
%Location=[a;b;c]--The coordinate of center at the bottom of borehole
%bore_dip=[direction;angle]--dip direction and angle of borehole
%H--deepness of borehole
%R--radius
%m--number of points on the circle
%
%h=figure;
t=linspace(0,2*pi,m);
x1=R*cos(t);
y1=R*sin(t);
z1=zeros(1,m);
z2=H*ones(1,m);
%coordinate transformation
alpha=pi/2-bore_dip(1,1);
phi=bore_dip(2,1);  
%transportation matrix, transport the local coordinate system into global coordinate system
nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
     sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
     -sin(phi),           0,           cos(phi);];
 for i=1:m
     Coordinate_down(:,i)=nto*[x1(i);y1(i);z1(i)]+location;
     Coordinate_up(:,i)=nto*[x1(i);y1(i);z2(i)]+location;
 end
colort=[0,0,0];
for i=1:m
    x=[Coordinate_down(1,i),Coordinate_up(1,i)];
    y=[Coordinate_down(2,i),Coordinate_up(2,i)];
    z=[Coordinate_down(3,i),Coordinate_up(3,i)];
    line(x,y,z,'Color',colort)
    hold on
end

for i=1:m-1
    x=[Coordinate_down(1,i),Coordinate_down(1,i+1)];
    y=[Coordinate_down(2,i),Coordinate_down(2,i+1)];
    z=[Coordinate_down(3,i),Coordinate_down(3,i+1)];
    line(x,y,z,'Color',colort)
    hold on
end
x=[Coordinate_down(1,m),Coordinate_down(1,1)];
y=[Coordinate_down(2,m),Coordinate_down(2,1)];
z=[Coordinate_down(3,m),Coordinate_down(3,1)];
line(x,y,z,'Color',colort)
hold on

for i=1:m-1
    x=[Coordinate_up(1,i),Coordinate_up(1,i+1)];
    y=[Coordinate_up(2,i),Coordinate_up(2,i+1)];
    z=[Coordinate_up(3,i),Coordinate_up(3,i+1)];
    line(x,y,z,'Color',colort)
    hold on
end
x=[Coordinate_up(1,m),Coordinate_up(1,1)];
y=[Coordinate_up(2,m),Coordinate_up(2,1)];
z=[Coordinate_up(3,m),Coordinate_up(3,1)];
line(x,y,z,'Color',colort)
hold on
return

