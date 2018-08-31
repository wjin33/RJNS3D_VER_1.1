function [theta,len_c_tra,len_d_tra,len_t_tra,plane_p]=GenRMWFT(x0,y0,z0,dx,dy,dz,joint_c,j_dd_dip,hm,rmc,mwdy,mwdx,mw_dd_dip)
%(x0,y0,z0)--coordinate of the first vertex of cuboid
%(dx,dy,dz)--length,width,hight of the cuboid
%hm--1xn matrix, handle of each joint
%j_dd_dip--2xn matrix, trend and dip of each joint
%joint_c--3xn matrix, center of each joint
%p_cn1--3x1 matrix, the global coordinates of left lower point of measuring window
%p_cn2--3x1 matrix, the global coordinates of right upper point ofmeasuring window
%is_jmw_flag--1:intersect£¬0:non-intersect
%plane_p--6xn matrix, coordinates of the endpoint of trace in local measure window coordinate system
%rmc--3xn matrix, center of each measure window
%mwdy--1xn  half length of the vertical axis of measure window
%mwdx--1xn  half length of the horizontal axis of measure window
%len_tc_tra=[]; len_td_tra=[]; len_tct_tra=[];
%c_tra-6xn; d_tra-6xn; t_tra-6xn;
%ratio_h--the ratio horizontal axis over radius,  ratio_v--the ratio vertical  axis over radius,
%theta--the angle bewteen trace and horizontal edge of the measure window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%genarate a global meansure window and translate the coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%transform to local measure windom coordinate system
cs_o=rmc;
alpha=pi/2-mw_dd_dip(1,1);
beta=mw_dd_dip(2,1);
nto=[cos(alpha)*cos(beta),-sin(alpha),cos(alpha)*sin(beta);
     sin(alpha)*cos(beta),cos(alpha), sin(alpha)*sin(beta);
     -sin(beta),0,cos(beta);];
otn=[cos(alpha)*cos(beta),sin(alpha)*cos(beta),-sin(beta);
     -sin(alpha),cos(alpha),0;
     cos(alpha)*sin(beta),sin(alpha)*sin(beta),cos(beta);];
%angular point of measure windom in local coordinate system
np_cn1=[mwdx;-mwdy;0];
np_cn2=[-mwdx;mwdy;0];
%transform to global coordinate system
p_cn1=nto*[dx;-dy;0]+cs_o; p_cn1=PruneMartix(p_cn1,1e-5);
p_cn2=nto*[-dx;dy;0]+cs_o; p_cn2=PruneMartix(p_cn2,1e-5);
%the four angular point in local coordinate system
nrec_cn1=np_cn1;
nrec_cn2=[np_cn1(1,1);np_cn2(2,1);0];
nrec_cn3=np_cn2;
nrec_cn4=[np_cn2(1,1);np_cn1(2,1);0];
%transform the four angular point in local coordinate system 
%to global coordinate system, anticlockwise
orec_cn1=nto*nrec_cn1+cs_o; 
orec_cn2=nto*nrec_cn2+cs_o; 
orec_cn3=nto*nrec_cn3+cs_o; 
orec_cn4=nto*nrec_cn4+cs_o; 
%draw the measure windom in global system
line([orec_cn1(1,1),orec_cn2(1,1)],[orec_cn1(2,1),orec_cn2(2,1)],[orec_cn1(3,1),orec_cn2(3,1)],'Color','k','LineWidth',2);
line([orec_cn2(1,1),orec_cn3(1,1)],[orec_cn2(2,1),orec_cn3(2,1)],[orec_cn2(3,1),orec_cn3(3,1)],'Color','k','LineWidth',2);
line([orec_cn3(1,1),orec_cn4(1,1)],[orec_cn3(2,1),orec_cn4(2,1)],[orec_cn3(3,1),orec_cn4(3,1)],'Color','k','LineWidth',2);
line([orec_cn4(1,1),orec_cn1(1,1)],[orec_cn4(2,1),orec_cn1(2,1)],[orec_cn4(3,1),orec_cn1(3,1)],'Color','k','LineWidth',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%genarate a global meansure window and the joint traces 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmw2=figure;%generate a new figure, 3-dimension
set(0,'CurrentFigure',hmw2)
colort=[0,0,0];
%subsurface, generate vertexes from the first, anticlockwise
x=[x0,x0+dx,x0+dx,x0+dx,x0+dx,x0,   x0,   x0];
y=[y0,y0,   y0,   y0+dy,y0+dy,y0+dy,y0+dy,y0];
z=[z0,z0,   z0,   z0,   z0,   z0,   z0,   z0];
line(x,y,z,'Color',colort)
hold on

%topsurface, generate vertexes from the first, anticlockwise
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
DisAxis([0;0;0],3);
view(3);
axis equal;
hold on
%draw the measure windom in global system
line([orec_cn1(1,1),orec_cn2(1,1)],[orec_cn1(2,1),orec_cn2(2,1)],[orec_cn1(3,1),orec_cn2(3,1)],'Color','k','LineWidth',2);
line([orec_cn2(1,1),orec_cn3(1,1)],[orec_cn2(2,1),orec_cn3(2,1)],[orec_cn2(3,1),orec_cn3(3,1)],'Color','k','LineWidth',2);
line([orec_cn3(1,1),orec_cn4(1,1)],[orec_cn3(2,1),orec_cn4(2,1)],[orec_cn3(3,1),orec_cn4(3,1)],'Color','k','LineWidth',2);
line([orec_cn4(1,1),orec_cn1(1,1)],[orec_cn4(2,1),orec_cn1(2,1)],[orec_cn4(3,1),orec_cn1(3,1)],'Color','k','LineWidth',2);
[row,col]=size(hm);
plane_p=[];
for i=1:col
    [np1,np2,is_jmw_flag]=IsJointMWPInsect(joint_c(:,i),j_dd_dip(:,i),hm(i),p_cn1,p_cn2,mw_dd_dip);
	if is_jmw_flag==1
		plane_p=[plane_p,[np1;np2];];   %store local coordinates of the intersected points
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%genarate a local meansure window and the joint traces 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hmw3=figure;%generate another new figure, 2-dimension
set(0,'CurrentFigure',hmw3)
%rotate 90 degree around z-axis
plane_p=[plane_p(2,:);-plane_p(1,:);plane_p(3,:);plane_p(5,:);-plane_p(4,:);plane_p(6,:)];
plane_p=PruneMartix(plane_p,1e-5);
%the four angular point in local coordinate system,which nrec_cn1 is original
nrec_cn1=nrec_cn1+[-dx;dy;0];
nrec_cn2=nrec_cn2+[-dx;dy;0];
nrec_cn3=nrec_cn3+[-dx;dy;0];
nrec_cn4=nrec_cn4+[-dx;dy;0];
%rotate 90 degree around z-axis
nrec_cn1=[nrec_cn1(2,:);-nrec_cn1(1,:);nrec_cn1(3,:)];
nrec_cn2=[nrec_cn2(2,:);-nrec_cn2(1,:);nrec_cn2(3,:)];
nrec_cn3=[nrec_cn3(2,:);-nrec_cn3(1,:);nrec_cn3(3,:)];
nrec_cn4=[nrec_cn4(2,:);-nrec_cn4(1,:);nrec_cn4(3,:)];
%draw the measure windom in local system
line([nrec_cn1(1,1),nrec_cn2(1,1)],[nrec_cn1(2,1),nrec_cn2(2,1)],[nrec_cn1(3,1),nrec_cn2(3,1)],'Color','k','LineWidth',2);
line([nrec_cn2(1,1),nrec_cn3(1,1)],[nrec_cn2(2,1),nrec_cn3(2,1)],[nrec_cn2(3,1),nrec_cn3(3,1)],'Color','k','LineWidth',2);
line([nrec_cn3(1,1),nrec_cn4(1,1)],[nrec_cn3(2,1),nrec_cn4(2,1)],[nrec_cn3(3,1),nrec_cn4(3,1)],'Color','k','LineWidth',2);
line([nrec_cn4(1,1),nrec_cn1(1,1)],[nrec_cn4(2,1),nrec_cn1(2,1)],[nrec_cn4(3,1),nrec_cn1(3,1)],'Color','k','LineWidth',2);
nvm=[nrec_cn1,nrec_cn2,nrec_cn3,nrec_cn4,nrec_cn1];
%set c_ inside the MW;  d_ divided by MW;  t_ cross the MW
c_tra=[]; d_tra=[]; t_tra=[];o_tra=[];
len_c_tra=[]; len_d_tra=[]; len_t_tra=[];
[row,col]=size(plane_p);
for i=1:col
    %lp_flag--0:inside the polygon, 1:line segment is cut, 2:line segment cross the ploygon, -1:outside the polygon
%rec_p--6x1 matrix -coordinates of the the coordinates of  intersection point and the pints inside the polygon

    [rec_p,lp_flag]=GetPolyLSIP2D(plane_p(1:3,i),plane_p(4:6,i),nvm);
	temp=plane_p(:,i);
	lens=sqrt(((temp(1,1)-temp(4,1))^2+(temp(2,1)-temp(5,1))^2));
	if lp_flag==0
		c_tra=[c_tra,plane_p(:,i)];
		len_c_tra=[len_c_tra,lens];
	elseif lp_flag==1
		d_tra=[d_tra,plane_p(:,i)];
		len_d_tra=[len_d_tra,lens];
	elseif lp_flag==2
		t_tra=[t_tra,plane_p(:,i)];
		len_t_tra=[len_t_tra,lens];
    elseif lp_flag==-1
        o_tra=[o_tra,plane_p(:,i)];
    end
end

[row,col]=size(t_tra);
if col~=0
  theta=atan(abs((t_tra(5,1)-t_tra(2,1))/(t_tra(4,1)-t_tra(1,1))));
    colort=[0,0,0]; 
	for i=1:col
		line(([t_tra(1,i),t_tra(4,i)]),([t_tra(2,i),t_tra(5,i)]),([t_tra(3,i),t_tra(6,i)]),'Color',colort,'LineWidth',2);
    end
end

[row,col]=size(d_tra);
if col~=0
    theta=atan(abs((d_tra(5,1)-d_tra(2,1))/(d_tra(4,1)-d_tra(1,1))));
    colort=[0,0.8,0];
	for i=1:col
		line(([d_tra(1,i),d_tra(4,i)]),([d_tra(2,i),d_tra(5,i)]),([d_tra(3,i),d_tra(6,i)]),'Color',colort,'LineWidth',2);
    end
end

[row,col]=size(c_tra);
if col~=0
     theta=atan(abs((c_tra(5,1)-c_tra(2,1))/(c_tra(4,1)-c_tra(1,1))));		
    colort=[1,0,0];
	for i=1:col
		line(([c_tra(1,i),c_tra(4,i)]),([c_tra(2,i),c_tra(5,i)]),([c_tra(3,i),c_tra(6,i)]),'Color',colort,'LineWidth',2);
    end
end

[row,col]=size(o_tra);
if col~=0
     theta=atan(abs((o_tra(5,1)-o_tra(2,1))/(o_tra(4,1)-o_tra(1,1))));
    colort=[0,0,0.8]; 
	for i=1:col
		line(([o_tra(1,i),o_tra(4,i)]),([o_tra(2,i),o_tra(5,i)]),([o_tra(3,i),o_tra(6,i)]),'Color',colort,'LineWidth',2);
    end
end
axis equal
end



