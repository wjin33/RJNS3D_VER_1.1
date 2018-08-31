function [mean_size,sigma]=CalmeanvarSize(lumtav,j_dd_dip,k_axis,gamma,bore_dip,radius,H,N_intersect,N_fullintsect)
%计算椭圆裂隙大小的均值和方差
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
b_ax_n1=[b_ax_o'*j_x_o; b_ax_o'*j_y_o; b_ax_o'*j_z_o;];%[I;J;K;]
miu=sqrt(b_ax_n1(2,1)^2+b_ax_n1(3,1)^2);
A=1/miu^2;
B=2*b_ax_n1(1,1)*b_ax_n1(2,1)/(b_ax_n1(3,1)*miu^2);
C=b_ax_n1(1,1)^2*b_ax_n1(2,1)^2/(b_ax_n1(3,1)^2*miu^2)+miu^2*k_axis^2/b_ax_n1(3,1)^2;
xi_a=sqrt(2/(A+C+sqrt((A-C)^2+B^2)));
xi_b=sqrt(2/(A+C-sqrt((A-C)^2+B^2)));
miu=sqrt(b_ax_n1(2,1)^2+b_ax_n1(3,1)^2);
denominator=sqrt((b_ax_n1(3,1)^2+b_ax_n1(1,1)^2*b_ax_n1(2,1)^2)/(k_axis^2*miu^2)+miu^2);
mean_size=(N_intersect-N_fullintsect)/(2*sqrt(2)*pi*lumtav*H*denominator);
sigma=sqrt((N_intersect/(pi*lumtav*H)-radius^2-sqrt(2)*mean_size*radius*denominator)/(xi_a*xi_b)-mean_size^2);
return
end











