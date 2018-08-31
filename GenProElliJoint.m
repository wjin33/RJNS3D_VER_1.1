function [h,coordinate]=GenProElliJoint(point_c,j_dd_dip,major_axis,k_axis,gamma,location,bore_dip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate the major and minor axis of the projected fracture in the
%%projection plane as well as the omega of retation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calcualte unit vector of borehole axis for coordinate system-1  in global coordinate system%%
alpha=pi/2-bore_dip(1,1);
phi=bore_dip(2,1);
%transportation matrix, transport the local coordinate system into global coordinate system
nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
     sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
     -sin(phi),           0,           cos(phi);];
b_ax_o=nto*[0.0;0.0;1.0]; %unit normal vector of the projection plane in global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=pi/2-j_dd_dip(1,1);  %the rotato angle, trend of joint substract trend of x-axis
beta=j_dd_dip(2,1);  %the dip angle
%transportation matrix, transport the local coordinate system-1 into global coordinate system
otn=[ cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma),  sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma), -sin(beta)*cos(gamma);
     -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma), -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma),  sin(beta)*sin(gamma);
      cos(alpha)*sin(beta),                                   sin(alpha)*sin(beta),                                   cos(beta);           ];
nto=otn';
j_x_o=nto*[1.0;0.0;0.0]; %unit x-axis vector of the joint in global coordinate system
j_y_o=nto*[0.0;1.0;0.0]; %unit y-axis vector of the joint plane in global coordinate system
j_z_o=nto*[0.0;0.0;1.0]; %unit normal vector of the joint plane in global coordinate system
b_ax_n1=[b_ax_o'*j_x_o; b_ax_o'*j_y_o; b_ax_o'*j_z_o;];
%transportation matrix, transport the local coordinate system-2 into the local coordinate system-1
miu=sqrt(b_ax_n1(2,1)^2+b_ax_n1(3,1)^2);
n1_t_n2=[miu,          -b_ax_n1(1,1)*b_ax_n1(2,1)/miu,   -b_ax_n1(1,1)*b_ax_n1(3,1)/miu;
         0,             b_ax_n1(3,1)/miu,                -b_ax_n1(2,1)/miu;
         b_ax_n1(1,1),  b_ax_n1(2,1),                     b_ax_n1(3,1);];
n2_t_n1=n1_t_n2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the projected center of the elllipse
location_n1=otn*(location-point_c);
n1_lc=[b_ax_n1(1,1);
       b_ax_n1(2,1);
       b_ax_n1(3,1);
       -(b_ax_n1')*location_n1;];  
lambada=-b_ax_o'*b_ax_o;
ppoint_c_n1=[b_ax_n1(1,1)*(n1_lc(4,1))/lambada;
             b_ax_n1(2,1)*(n1_lc(4,1))/lambada;
             b_ax_n1(3,1)*(n1_lc(4,1))/lambada;];
ppoint_c=nto*ppoint_c_n1+point_c;
          
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Drawing the local coordinate system-1
% cs_o=point_c;
% x=[cs_o(1,1),cs_o(1,1)+25*j_x_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+25*j_x_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+25*j_x_o(3,1)];
% line(x,y,z,'Color','k')
% hold on
% x=[cs_o(1,1),cs_o(1,1)+5*j_y_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+5*j_y_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+5*j_y_o(3,1)];
% line(x,y,z,'Color','g')
% hold on
% x=[cs_o(1,1),cs_o(1,1)+5*j_z_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+5*j_z_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+5*j_z_o(3,1)];
% line(x,y,z,'Color','r')
% hold on
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Drawing the local coordinate system-2
% pj_x_o=nto*n2_t_n1*[1.0;0.0;0.0]; %unit x-axis vector of the joint in global coordinate system
% pj_y_o=nto*n2_t_n1*[0.0;1.0;0.0]; %unit y-axis vector of the joint plane in global coordinate system
% pj_z_o=nto*n2_t_n1*[0.0;0.0;1.0]; %unit normal vector of the joint plane in global coordinate system
% cs_o=ppoint_c;
% x=[cs_o(1,1),cs_o(1,1)+25*pj_x_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+25*pj_x_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+25*pj_x_o(3,1)];
% line(x,y,z,'Color','b')
% hold on
% x=[cs_o(1,1),cs_o(1,1)+5*pj_y_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+5*pj_y_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+5*pj_y_o(3,1)];
% line(x,y,z,'Color','b')
% hold on
% x=[cs_o(1,1),cs_o(1,1)+5*pj_z_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+5*pj_z_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+5*pj_z_o(3,1)];
% line(x,y,z,'Color','b')
% hold on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=1/miu^2;
B=2*b_ax_n1(1,1)*b_ax_n1(2,1)/(b_ax_n1(3,1)*miu^2);
C=b_ax_n1(1,1)^2*b_ax_n1(2,1)^2/(b_ax_n1(3,1)^2*miu^2)+miu^2*k_axis^2/b_ax_n1(3,1)^2;
xi_a=sqrt(2/(A+C-sqrt((A-C)^2+B^2)));
xi_b=sqrt(2/(A+C+sqrt((A-C)^2+B^2)));
pmajor_axis=xi_a*major_axis/2;
pminor_axis=xi_b*major_axis/2;
if A==C
    omega=pi/4;
else
    omega=atan(B/(A-C))/2;
end
sinomega = sin(omega);
cosomega = cos(omega);
angle = linspace(0, 360, 40).* (pi / 180);
sinalpha = sin(angle);
cosalpha = cos(angle);
[~,col]=size(angle);
pj_v_n2=zeros(3,col);
for i=1:col
    pj_v_n2(1,i) = pmajor_axis * cosalpha(1,i) * cosomega - pminor_axis * sinalpha(1,i) * sinomega;
    pj_v_n2(2,i) = pmajor_axis * cosalpha(1,i) * sinomega + pminor_axis * sinalpha(1,i) * cosomega;
end
pj_v_n2=PruneMartix(pj_v_n2,1e-5);

pj_v_o=zeros(3,col);
[~,col]=size(pj_v_n2);
t_m=nto*n2_t_n1*pj_v_n2;
for m=1:col
	pj_v_o(:,m)=t_m(:,m)+ppoint_c;
end
colort=rand(1,3);
h=patch('Vertices',(pj_v_o'),'Faces',(1:1:col),'FaceColor',colort,'EdgeColor',colort);
coordinate=pj_v_o;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following are some failed codes, but they are useful in some latter
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%求两个椭圆的相交函数
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cs_o=point_c;
% alpha=pi/2-j_dd_dip(1,1);  %the rotato angle, trend of joint substract trend of x-axis
% phi=j_dd_dip(2,1);  %the rotato angle of dip 
% %transportation matrix, transport the local coordinate system into global coordinate system
% nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
%      sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
%      -sin(phi),           0,           cos(phi);];
% j_ov_n=nto*[0.0;0.0;10]; % normal vector of the joint plane in global coordinate system
% j_ov_x=nto*[1.0;0.0;0.0];% x-axis vector of the joint plane in global coordinate system
% j_ov_y=nto*[0.0;1.0;0.0];% y-axis vector of the joint plane in global coordinate system
% o_lc=[j_ov_n(1,1);j_ov_n(2,1);j_ov_n(3,1);-(j_ov_n')*cs_o;];  %generate coefficients vector of the joint plane equation (oax+oby+ocz+od=0) in global coordinate system
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%  point_o_bore=Location;
%  point_o_bore(3,1)=-(Location(1,1)*o_lc(1,1)+Location(2,1)*o_lc(2,1)+o_lc(4,1))/o_lc(3,1);
%  p_n_x=(cs_o-point_o_bore)'*j_ov_x;
% %  p_o_interceptx=[point_o_bore(1,1)+Lambada*j_ov_x(1,1);
% %                  point_o_bore(2,1)+Lambada*j_ov_x(2,1);
% %                  point_o_bore(3,1)+Lambada*j_ov_x(3,1)];
%  p_n_y=(cs_o-point_o_bore)'*j_ov_y;
% %  p_o_intercepty=[point_o_bore(1,1)+Lambada*j_ov_y(1,1);
% %                  point_o_bore(2,1)+Lambada*j_ov_y(2,1);
% %                  point_o_bore(3,1)+Lambada*j_ov_y(3,1)];
% 
% 
% theta=0:0.005:pi;
% %parametric equation for an ellipse which is formed by interception between
% %borehole wall and the joint plane
% bore_ellipse=[Radius_c*cos(theta)/cos(phi);Radius_c*sin(theta)];
% %parametric equation for the joint in local coordinate system-2
% joint_ellipse=[-p_n_x+major_axis*cos(gamma)*cos(theta)+major_axis/(k_axis)*sin(gamma)*sin(theta);
%                 p_n_y-major_axis*sin(gamma)*cos(theta)+major_axis/(k_axis)*cos(gamma)*sin(theta)];
% [row,col]=size(theta);
% for i=1:col
%         e_nvm(1,i) = major_axis/2 * cosalpha(1,i) * cosbeta - major_axis/(2*k_axis) * sinalpha(1,i) * sinbeta;
%         e_nvm(2,i) = major_axis/2 * cosalpha(1,i) * sinbeta + major_axis/(2*k_axis) * sinalpha(1,i) * cosbeta;
% end

     
 