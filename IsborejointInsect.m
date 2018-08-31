function [flag_int,flag_full]=IsborejointInsect(point_c,j_dd_dip,major_axis,k_axis,gamma,location,bore_dip,radius)
%%ÅÐ¶ÏÁÑÏ¶Óë×ê¿×½»ÇÐÇé¿ö
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
pro_bo_point_c=location;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=1/miu^2;
B=2*b_ax_n1(1,1)*b_ax_n1(2,1)/(b_ax_n1(3,1)*miu^2);
C=b_ax_n1(1,1)^2*b_ax_n1(2,1)^2/(b_ax_n1(3,1)^2*miu^2)+miu^2*k_axis^2/b_ax_n1(3,1)^2;
xi_a=sqrt(2/(A+C+sqrt((A-C)^2+B^2)));
xi_b=sqrt(2/(A+C-sqrt((A-C)^2+B^2)));
if xi_a>=xi_b
    pmajor_axis=xi_a*major_axis/2;
    pminor_axis=xi_b*major_axis/2;
else
    pmajor_axis=xi_b*major_axis/2;
    pminor_axis=xi_a*major_axis/2;
end
if A==C
    omega=pi/4;
else
    omega=atan(B/(A-C))/2;
end
n3_t_n2=[cos(omega),   -sin(omega),      0;
         sin(omega),    cos(omega),      0;
         0,             0,               1;];

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Drawing the local coordinate system-3
% pj_x_o=nto*n2_t_n1*n3_t_n2*[1.0;0.0;0.0]; %unit x-axis vector of the joint in global coordinate system
% pj_y_o=nto*n2_t_n1*n3_t_n2*[0.0;1.0;0.0]; %unit y-axis vector of the joint plane in global coordinate system
% pj_z_o=nto*n2_t_n1*n3_t_n2*[0.0;0.0;1.0]; %unit normal vector of the joint plane in global coordinate system
% cs_o=pro_bo_point_c;
% %cs_o=ppoint_c;
% x=[cs_o(1,1),cs_o(1,1)+25*pj_x_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+25*pj_x_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+25*pj_x_o(3,1)];
% line(x,y,z,'Color','r')
% hold on
% x=[cs_o(1,1),cs_o(1,1)+5*pj_y_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+5*pj_y_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+5*pj_y_o(3,1)];
% line(x,y,z,'Color','r')
% hold on
% x=[cs_o(1,1),cs_o(1,1)+5*pj_z_o(1,1)];
% y=[cs_o(2,1),cs_o(2,1)+5*pj_z_o(2,1)];
% z=[cs_o(3,1),cs_o(3,1)+5*pj_z_o(3,1)];
% line(x,y,z,'Color','r')
% hold on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Drawing the specifical curve_L
% sinomega = sin(omega);
% cosomega = cos(omega);
% angle = linspace(0, 360, 50).* (pi / 180);
% sinalpha = sin(angle);
% cosalpha = cos(angle);
% [~,col]=size(angle);
% pj_v_n2=zeros(3,col);
% for i=1:col
%     pj_v_n2(1,i) = pmajor_axis * cosalpha(1,i) * cosomega - pminor_axis * sinalpha(1,i) * sinomega +  ...
%                   (pminor_axis * cosalpha(1,i) * cosomega - pmajor_axis * sinalpha(1,i) * sinomega)*radius/sqrt(pminor_axis^2*cosalpha(1,i).^2+pmajor_axis^2*sinalpha(1,i).^2);
%     pj_v_n2(2,i) = pmajor_axis * cosalpha(1,i) * sinomega + pminor_axis * sinalpha(1,i) * cosomega + ...
%                   (pminor_axis * cosalpha(1,i) * sinomega + pmajor_axis * sinalpha(1,i) * cosomega)*radius/sqrt(pminor_axis^2*cosalpha(1,i).^2+pmajor_axis^2*sinalpha(1,i).^2);
% end
% pj_v_n2=PruneMartix(pj_v_n2,1e-5);
% 
% pj_v_o=zeros(3,col);
% [~,col]=size(pj_v_n2);
% t_m=nto*n2_t_n1*pj_v_n2;
% for m=1:col
% 	pj_v_o(:,m)=t_m(:,m)+pro_bo_point_c;
% end
% colort=rand(1,3);
% h=patch('Vertices',(pj_v_o'),'Faces',(1:1:col),'FaceColor',colort,'EdgeColor',colort);

%%

angle = linspace(-180, 180, 360).* (pi / 180);
sinalpha = sin(angle);
cosalpha = cos(angle);
[~,col]=size(angle);
pj_v_n3=zeros(3,col);
for i=1:col
    pj_v_n3(1,i) = pmajor_axis * cosalpha(1,i) + (pminor_axis * radius * cosalpha(1,i))/sqrt(pminor_axis^2*cosalpha(1,i).^2+pmajor_axis^2*sinalpha(1,i).^2);
    pj_v_n3(2,i) = pminor_axis * sinalpha(1,i) + (pmajor_axis * radius * sinalpha(1,i))/sqrt(pminor_axis^2*cosalpha(1,i).^2+pmajor_axis^2*sinalpha(1,i).^2);
end
pj_v_n3=PruneMartix(pj_v_n3,1e-5);
diff_theta=zeros(1,col);
theta2=zeros(1,col);
for i=1:col
    if i<=col/2
        theta2(1,i)=-acos([1,0,0]*pj_v_n3(:,i)/sqrt(pj_v_n3(1,i)^2+pj_v_n3(2,i)^2));
    else
        theta2(1,i)=acos([1,0,0]*pj_v_n3(:,i)/sqrt(pj_v_n3(1,i)^2+pj_v_n3(2,i)^2));
    end
    diff_theta(1,i)=angle(1,i)-theta2(1,i);
end
cen_bj_o=ppoint_c-pro_bo_point_c;   %vector of the line formed bewteen projected borehole center and projected fracture center in global coordinate system
dis_bcen_jcen=sqrt(cen_bj_o(1,1)^2+cen_bj_o(2,1)^2+cen_bj_o(3,1)^2);
pj_x_o=nto*n2_t_n1*n3_t_n2*[1.0;0.0;0.0];   
zeta=acos(cen_bj_o'*pj_x_o/dis_bcen_jcen);
pj_y_o=nto*n2_t_n1*n3_t_n2*[0.0;1.0;0.0];  
eta=acos(cen_bj_o'*pj_y_o/dis_bcen_jcen);
if eta>=pi/2;
   zeta=-zeta;
end
for i=1:col-1
    if zeta>=theta2(1,i) && zeta<theta2(1,i+1)
        eta=zeta+diff_theta(1,i);
        break
    end
end
pj_v_n3(1,1) = pmajor_axis * cos(eta) + (pminor_axis * radius * cos(eta) )/sqrt(pminor_axis^2*cos(eta).^2+pmajor_axis^2*sin(eta).^2);
pj_v_n3(2,1) = pminor_axis * sin(eta) + (pmajor_axis * radius * sin(eta) )/sqrt(pminor_axis^2*cos(eta).^2+pmajor_axis^2*sin(eta).^2);

%is_point=nto*n2_t_n1*n3_t_n2*[pj_v_n3(1,1);pj_v_n3(2,1);0.0]+pro_bo_point_c
% x=[pro_bo_point_c(1,1),ppoint_c(1,1)];
% y=[pro_bo_point_c(2,1),ppoint_c(2,1)];
% z=[pro_bo_point_c(3,1),ppoint_c(3,1)];
% line(x,y,z,'color','r','line','-')
% x=[pro_bo_point_c(1,1),is_point(1,1)];
% y=[pro_bo_point_c(2,1),is_point(2,1)];
% z=[pro_bo_point_c(3,1),is_point(3,1)];
% line(x,y,z,'color','k','line','--')

dis_bcen_jten=sqrt(pj_v_n3(1,1)^2+pj_v_n3(2,1)^2);
if dis_bcen_jcen<=dis_bcen_jten
    flag_int=1;
else
    flag_int=0;
end

%%

% angle = linspace(-180, 180, 360).* (pi / 180);
% sinalpha = sin(angle);
% cosalpha = cos(angle);
% [~,col]=size(angle);
% pj_v_n3=zeros(3,col);
for i=1:col
    pj_v_n3(1,i) = pmajor_axis * cosalpha(1,i) - (pminor_axis * radius * cosalpha(1,i))/sqrt(pminor_axis^2*cosalpha(1,i).^2+pmajor_axis^2*sinalpha(1,i).^2);
    pj_v_n3(2,i) = pminor_axis * sinalpha(1,i) - (pmajor_axis * radius * sinalpha(1,i))/sqrt(pminor_axis^2*cosalpha(1,i).^2+pmajor_axis^2*sinalpha(1,i).^2);
end
pj_v_n3=PruneMartix(pj_v_n3,1e-5);
diff_theta=zeros(1,col);
theta2=zeros(1,col);
for i=1:col
    if i<=col/2
        theta2(1,i)=-acos([1,0,0]*pj_v_n3(:,i)/sqrt(pj_v_n3(1,i)^2+pj_v_n3(2,i)^2));
    else
        theta2(1,i)=acos([1,0,0]*pj_v_n3(:,i)/sqrt(pj_v_n3(1,i)^2+pj_v_n3(2,i)^2));
    end
    diff_theta(1,i)=angle(1,i)-theta2(1,i);
end
% cen_bj_o=ppoint_c-pro_bo_point_c;   %vector of the line formed bewteen projected borehole center and projected fracture center in global coordinate system
% dis_bcen_jcen=sqrt(cen_bj_o(1,1)^2+cen_bj_o(2,1)^2+cen_bj_o(3,1)^2);
% pj_x_o=nto*n2_t_n1*n3_t_n2*[1.0;0.0;0.0];   
% zeta=acos(cen_bj_o'*pj_x_o/dis_bcen_jcen);
% pj_y_o=nto*n2_t_n1*n3_t_n2*[0.0;1.0;0.0];  
% eta=acos(cen_bj_o'*pj_y_o/dis_bcen_jcen);
% if eta>=pi/2;
%    zeta=-zeta;
% end
for i=1:col-1
    if zeta>=theta2(1,i) && zeta<theta2(1,i+1)
        eta=zeta+diff_theta(1,i);
        break
    end
end
pj_v_n3(1,1) = pmajor_axis * cos(eta) - (pminor_axis * radius * cos(eta) )/sqrt(pminor_axis^2*cos(eta).^2+pmajor_axis^2*sin(eta).^2);
pj_v_n3(2,1) = pminor_axis * sin(eta) - (pmajor_axis * radius * sin(eta) )/sqrt(pminor_axis^2*cos(eta).^2+pmajor_axis^2*sin(eta).^2);
dis_bcen_jten=sqrt(pj_v_n3(1,1)^2+pj_v_n3(2,1)^2);
if dis_bcen_jcen<=dis_bcen_jten
    flag_full=1;
else
    flag_full=0;
end
return


