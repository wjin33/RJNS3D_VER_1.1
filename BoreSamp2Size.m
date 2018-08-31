clear all;close all;
x0=0;y0=0;z0=0;
dx=80;dy=80;dz=80;
%lumta--density
lumtav=0.05;
point_c=PoissonProcess3D2(lumtav,dx,dy,dz);
%point_c=[dx/4;dy/4;dz/2];
[row,num_joint]=size(point_c);
%Uniformly distributed characteristic dimension   B1=0.5;B2=8;
major_axis = exprnd(20,1,num_joint);  %长轴而不是长半轴
% major_axis=unifrnd(10,60,1,num_joint);
% major_axis=30*ones(1,num_joint);


%One joint set with a deterministic orientation theta=pi/3; phi=pi/3;
j_dd_dip=[(0/3.0)*pi*ones(1,num_joint);(1.0/6.0)*pi*ones(1,num_joint)];
%generate the simulation region
[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,dx,dy,dz);
DisAxis([0;0;0;],3); % the global coordinate red--z; blue--y; green--x;
view(3);  axis equal
patch_handle=zeros(2,num_joint);
% the determined parameters of elliptical joints
gamma=0*pi/6;
k_axis=2;
% generate each elliptical joint inside the simulation region one by one

% for i=1:num_joint
%     patch_handle(1,i)=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,gamma,r_v_m);
% end

%number of boreholes drilled
num_borehole=5;
% bore hole surface location
location=[40,40,40,40,40;
          20,30,40,50,60;
          0,0,0,0,0;];
%borehole axis direction--vector
bore_dip=[1*pi/4*ones(1,num_borehole);0*pi/4*ones(1,num_borehole);];
%borehole deepness
H=80*ones(1,num_borehole);
%m-圆周方向的等分数
m=20;
%Radii
radius=1*ones(1,num_borehole);
for j=1:num_borehole
    Genborehole(location(:,j),bore_dip(:,j),H(:,j),radius(:,j),m);
end

flag_int=zeros(num_joint,num_borehole);
flag_full=zeros(num_joint,num_borehole);
for i=1:num_joint
    for j=1:num_borehole
%         patch_handle(2,i)=GenProElliJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,gamma,location(:,j),bore_dip(:,j));
        [flag_int(i,j),flag_full(i,j)]=IsborejointInsect(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,gamma,location(:,j),bore_dip(:,j),radius(:,j));
    end
end
N_intersect=sum(flag_int);
N_fullintsect=sum(flag_full);
mean_size=zeros(num_borehole,1);
sigma=zeros(num_borehole,1);
for i=1:num_borehole
    [mean_size(i,1),sigma(i,1)]=CalmeanvarSize(lumtav,j_dd_dip(:,1),k_axis,gamma,bore_dip(:,i),radius(:,i),H(:,i),N_intersect(1,i),N_fullintsect(1,i));
end
% 
sigma=mean(sigma);
mean_size=mean(mean_size);

% 
t_intermediate=sum(flag_int,2);
N1=0;N2=0;N3=0;N4=0;N5=0;
for i=1:num_joint
    if t_intermediate(i,1)==1
        N1=N1+1;
    elseif t_intermediate(i,1)==2
        N2=N2+1;
    elseif t_intermediate(i,1)==3
        N3=N3+1;
    elseif t_intermediate(i,1)==4
        N4=N4+1;
    elseif t_intermediate(i,1)==5
        N5=N5+1;
    end
end


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
A=1/miu^2;
B=2*b_ax_n1(1,1)*b_ax_n1(2,1)/(b_ax_n1(3,1)*miu^2);
C=b_ax_n1(1,1)^2*b_ax_n1(2,1)^2/(b_ax_n1(3,1)^2*miu^2)+miu^2*k_axis^2/b_ax_n1(3,1)^2;
xi_a=sqrt(2/(A+C-sqrt((A-C)^2+B^2)));
xi_b=sqrt(2/(A+C+sqrt((A-C)^2+B^2)));

save data_mul N1 N2 N3 N4 N5

%%

D_a=(2*mean_size-a_min)/(mean_size-a_min);
radius=radius(1,1);
a_min=xi_a*radius/xi_b^2;
D=10;%Distance between two boreholes
range=zeros(num_borehole,1);
for i=1:num_borehole
    range(i,1)=(i*D/2-radius)/xi_a;
end

%%

N_cal=zeros(num_borehole,1);
%a=xi_a*radius/xi_b^2:0.001:range(1,1);
a=a_min:0.05:range(1,1);
denominator=sqrt((b_ax_n1(3,1)^2+b_ax_n1(1,1)^2*b_ax_n1(2,1)^2)/(k_axis^2*miu^2)+miu^2);
fun_S1=(pi*xi_a*xi_b.*a.^2+pi*radius^2+pi*sqrt(2)*denominator.*a);
% fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S=num_borehole*fun_S1;
integral_1=trapz(a,fun_a.*fun_S);

i=2;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S1=(pi*xi_a*xi_b.*a.^2+pi*radius^2+pi*sqrt(2)*denominator.*a);
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=num_borehole*fun_S1-2*(num_borehole-1)*fun_S2;
integral_2=trapz(a,fun_a.*fun_S);

i=3;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S1=(pi*xi_a*xi_b.*a.^2+pi*radius^2+pi*sqrt(2)*denominator.*a);
fun_S=num_borehole*fun_S1-2*(num_borehole-1)*fun_S2+(num_borehole-2)*fun_S3;
integral_3=trapz(a,fun_a.*fun_S);  

i=4;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
% fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S1=(pi*xi_a*xi_b.*a.^2+pi*radius^2+pi*sqrt(2)*denominator.*a);
fun_S=num_borehole*fun_S1-2*(num_borehole-1)*fun_S2+(num_borehole-2)*fun_S3;
integral_4=trapz(a,fun_a.*fun_S);

i=5;
a=range(i-1,1):0.1:100;
%   fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
% fun_S5=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
% fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S1=(pi*xi_a*xi_b.*a.^2+pi*radius^2+pi*sqrt(2)*denominator.*a);

fun_S=num_borehole*fun_S1-2*(num_borehole-1)*fun_S2+(num_borehole-2)*fun_S3;
integral_5=trapz(a,fun_a.*fun_S);

N_cal(1,1)=lumtav*H(1,1)*(integral_1+integral_2++integral_3+integral_4+integral_5);




%%

i=2;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-1)*fun_S2;
integral_2=trapz(a,fun_a.*fun_S);
i=3;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-1)*fun_S2-2*(num_borehole-2)*fun_S3;
integral_3=trapz(a,fun_a.*fun_S);  
i=4;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-1)*fun_S2-2*(num_borehole-2)*fun_S3+(num_borehole-3)*fun_S4;
integral_4=trapz(a,fun_a.*fun_S);
i=5;
a=range(i-1,1):0.1:100;
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
% fun_S5=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S2=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-1)*fun_S2-2*(num_borehole-2)*fun_S3+(num_borehole-3)*fun_S4;%+(num_borehole-4)*fun_S5;
integral_5=trapz(a,fun_a.*fun_S);

N_cal(2,1)=lumtav*H(1,1)*(integral_2+integral_3+integral_4+integral_5);



%%

i=3;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-2)*fun_S3;
integral_3=trapz(a,fun_a.*fun_S);  
i=4;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-2)*fun_S3-2*(num_borehole-3)*fun_S4;
integral_4=trapz(a,fun_a.*fun_S);
i=5;
a=range(i-1,1):0.1:100;
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S5=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S3=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-2)*fun_S3-2*(num_borehole-3)*fun_S4+(num_borehole-4)*fun_S5;
integral_5=trapz(a,fun_a.*fun_S);

N_cal(3,1)=lumtav*H(1,1)*(integral_3+integral_4+integral_5);



%%

i=4;
a=range(i-1,1):0.1:range(i,1);
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-3)*fun_S4;
integral_4=trapz(a,fun_a.*fun_S);
i=5;
a=range(i-1,1):0.1:100;
%  fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S5=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
i=i-1;
fun_S4=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
fun_S=(num_borehole-3)*fun_S4-2*(num_borehole-4)*fun_S5;
integral_5=trapz(a,fun_a.*fun_S);

N_cal(4,1)=lumtav*H(1,1)*(integral_4+integral_5);


%%
i=5;
a=range(i-1,1):0.1:100;
% fun_a=0.1*exp(-0.1.*a);
% fun_a=lognpdf(a,10,sigma);
fun_a=(D_a-1)*a_min^(D_a-1).*a.^-D_a;
fun_S5=2*(xi_a*xi_b.*a.^2+radius*(xi_a+xi_b).*a+radius^2).*(acos((i-1)*D./(2*xi_a.*a+2*radius))-(i-1)*D*sqrt(1-(i-1)^2*D^2./(4*(xi_a.*a+radius).^2))./(2*xi_a.*a+2*radius));
integral_5=trapz(a,fun_a.*fun_S5);

N_cal(5,1)=lumtav*H(1,1)*integral_5;


% save data_exp_mul N_cal
% save data_lognormal_mul N_cal
save data_fractal_mul N_cal

%%postprocessor
% figure
% load data_mul
% N=[N1,N2,N3,N4,N5];
% bar(N)
% hold on
% a=1:5;
% b=zeros(5,3);
% load data_exp_mul
% b(:,1)=N_cal;
% load data_lognormal_mul
% b(:,2)=N_cal;
% load data_fractal_mul
% b(:,3)=N_cal;
% plot(a,b(:,1),'o',a,b(:,2),'*',a,b(:,3),'s')
% legend('Actual','Exponential','Lognormal','Fractal')
        




















