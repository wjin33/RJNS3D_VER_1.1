clear
clc
%??????
x0=0;
y0=0;
z0=0;
dx=80;
dy=80;
dz=80;
%lumta?????
lumtav=0.0006;
point_c=PoissonProcess3D2(lumtav,dx,dy,dz);
[row,num_joint]=size(point_c);
%?????0-8?????,???4?E(D)=8,E(lumtav)=lumtav,E(lumtaa)=lumtav*E(D)
radius=4*ones(1,num_joint);
j_dd_dip=[zeros(1,num_joint);(1.0/3.0)*pi*ones(1,num_joint)];
[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,dx,dy,dz);
DisAxis([0;0;0;],3); view(3);  axis equal
patch_handle=zeros(1,num_joint);
for i=1:num_joint
patch_handle(1,i)=GenDiskJoint(point_c(:,i),j_dd_dip(:,i),radius(i),r_v_m);
end

%????
mwh=10;  mwv=10;
[r_v_m,hcuoid]=GenCuboid((0.5*(dx-2*mwh)),(0.5*(dy-2*mwv)),(0.5*(dy-2*mwv)),(2*mwh),(2*mwh),(2*mwh));
%DisAxis([0;0;0;],3); 
view(3);  axis equal
cpoint_c=[];
cj_dd_dip=[];
cradius=[];
cpatch_handle=[];
for i=1:num_joint
h=GenDiskJoint(point_c(:,i),j_dd_dip(:,i),radius(i),r_v_m);
%??????????????
if h~=-1
cpatch_handle=[cpatch_handle,h];
cpoint_c=[cpoint_c,point_c(:,i)];
cj_dd_dip=[cj_dd_dip,j_dd_dip(:,i)];
cradius=[cradius,radius(i)];
end
end
