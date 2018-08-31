clear all;close all;
x0=0;y0=0;z0=0;
dx=80;dy=80;dz=80;
%lumta--density
lumtav=0.00001;
%point_c=PoissonProcess3D2(lumtav,dx,dy,dz);
point_c=[dx/2;dy/2;dz/2];
[row,num_joint]=size(point_c);
%Uniformly distributed characteristic dimension   B1=0.5;B2=8;
%major_axis = exprnd(10,1,num_joint);  %长轴而不是长半轴
% major_axis=unifrnd(10,60,1,num_joint);
 major_axis=50*ones(1,num_joint);

%One joint set with a deterministic orientation theta=pi/3; phi=pi/3;
%j_dd_dip=[(51/180)*pi*ones(1,num_joint);(50/180)*pi*ones(1,num_joint)];
j_dd_dip=[(158/180)*pi*ones(1,num_joint);(42/180)*pi*ones(1,num_joint)];
%generate the simulation region
[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,dx,dy,dz);
DisAxis([0;0;0;],3); % the global coordinate red--z; blue--y; green--x;
view(3);  axis equal
patch_handle=zeros(2,num_joint);
% the determined parameters of elliptical joints
for angle=1:180
gamma=(angle/180)*pi;
k_axis=1.5;
% generate each elliptical joint inside the simulation region one by one

% for i=1:num_joint
%     patch_handle(1,i)=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,gamma,r_v_m);
% end

%number of boreholes drilled
num_borehole=5;
% bore hole surface location
location=[20,30,40,50,60;
          40,40,40,40,40;
          0,0,0,0,0;];
%borehole axis direction--vector
bore_dip=[1*pi/4*ones(1,num_borehole);0*pi/4*ones(1,num_borehole);];
%borehole deepness
H=80*ones(1,num_borehole);
%m-圆周方向的等分数
m=20;
%Radii
radius=3*ones(1,num_borehole);
for j=1:num_borehole
    Genborehole(location(:,j),bore_dip(:,j),H(:,j),radius(:,j),m);
end

for i=1:num_joint
    for j=1:num_borehole
         [patch_handle(2,i),coor_pro_elli]=GenProElliJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,gamma,location(:,j),bore_dip(:,j));
    end
end

com_para_x=(max(coor_pro_elli(1,:))+min(coor_pro_elli(1,:)))/2;
[com_para_y, indice]=max(coor_pro_elli(2,:));
if abs(com_para_x-coor_pro_elli(1,indice))<0.5 && max(coor_pro_elli(1,:))>=com_para_y
    break
end
 end

angle


