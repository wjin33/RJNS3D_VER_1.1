clear all; close all;
x0=0;y0=0;z0=0;
dx=10;dy=10;dz=10;
%lumta--density
lumtav=0.2;
%rotation angle
beta=0*pi/2;
%axis ratio
k_axis=1;
%fracture center
point_c=PoissonProcess3D2(lumtav,dx,dy,dz);
[row,num_joint]=size(point_c);
%Dip direction;Dip angle
j_dd_dip=pi/4*ones(num_joint,2);
%fracture orientation--Fisher distribution
% es_k=9;%Fisher parameter
% j_dd_dip=GenFisherRand(jg_dd_dip,es_k,num_joint);
%Normally distributed characteristic dimension  with mean=0.7939 and
%standard difference=1
% half_major_axis = normrnd(0.7939,1,1,num_joint);
major_axis=0.2*ones(1,num_joint);

%generate the simulation region
[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,dx,dy,dz);
DisAxis([0;0;0;],3); % the global coordinate red--z; blue--y; green--x;
view(3);  axis equal
patch_handle=zeros(1,num_joint);
% the determined parameters of elliptical joints


% generate each elliptical joint inside the simulation region one by one
for i=1:num_joint
patch_handle(1,i)=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,r_v_m);
end

break


%Measure window #1
rmc=[2.5*cos(1*pi/180);2.5*sin(1*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[359*pi/180.0;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc1,theta1,len_c_tra1,len_d_tra1,len_t_tra1,tra_coord1]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);


%Measure window #2
rmc=[2.5*cos(15*pi/180);2.5*sin(15*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[345*pi/180.0;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc2,theta2,len_c_tra2,len_d_tra2,len_t_tra2,tra_coord2]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);

%Measure window #3
rmc=[2.5*cos(30*pi/180);2.5*sin(30*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[330*pi/180;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc3,theta3,len_c_tra3,len_d_tra3,len_t_tra3,tra_coord3]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);

%Measure window #4
rmc=[2.5*cos(45*pi/180);2.5*sin(45*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[315*pi/180;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc4,theta4,len_c_tra4,len_d_tra4,len_t_tra4,tra_coord4]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);

%Measure window #5
rmc=[2.5*cos(60*pi/180);2.5*sin(60*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[300*pi/180;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc5,theta5,len_c_tra5,len_d_tra5,len_t_tra5,tra_coord5]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);

%Measure window #6
rmc=[2.5*cos(75*pi/180);2.5*sin(75*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[285*pi/180;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc6,theta6,len_c_tra6,len_d_tra6,len_t_tra6,tra_coord6]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);


%Measure window #7
rmc=[2.5*cos(89*pi/180);2.5*sin(89*pi/180);0.5*dz];% the center coordinates of measuring widow
mw_dd_dip=[271*pi/180;1*pi/2.0]; % the orientation of measuring widow
mwh=2.5;  mwv=5; % the vertical and horizontal length of measuring widow
[beta_loc7,theta7,len_c_tra7,len_d_tra7,len_t_tra7,tra_coord7]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);
