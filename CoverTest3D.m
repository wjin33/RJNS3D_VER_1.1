function cv_flag=CoverTest3D(sx0,sy0,sz0,si,sj,sk,sx,sy,sz,h,point_c,j_dd_dip,~)
%cv_flag--1:cover,  -1:non-cover
%sx0,sy0,sz0--initial point of the model
%si,sj,sk--the index of each box
%sx,sy,sz-scale in 3 directions
%h--handle of this joint
%point_c--3x1 matrix, center of the joint
%j_dd_dip--2x1 matrix, strike of the joint
o_j_p=get(h,'Vertices');
o_j_p=o_j_p';

x0=sx0+(si-1)*sx;  y0=sy0+(sj-1)*sy;  z0=sz0+(sk-1)*sz;
r_v_m=zeros(24,5);

%face 1, the coordinates and attitude of topsurface 
x=x0; y=y0; z=z0+sz;
r_v_m(1,1:3)=[x,y,z];
x=x0+sx; y=y0; z=z0+sz;
r_v_m(2,1:3)=[x,y,z];
x=x0+sx; y=y0+sy; z=z0+sz;
r_v_m(3,1:3)=[x,y,z];
x=x0; y=y0+sy; z=z0+sz;
r_v_m(4,1:3)=[x,y,z];
%dip and trend
for i=1:4
    r_v_m(i,4:5)=[pi/2.0,0];
end

%face 2
x=x0; y=y0+sy; z=z0+sz;
r_v_m(5,1:3)=[x,y,z];
x=x0; y=y0+sy; z=z0;
r_v_m(6,1:3)=[x,y,z];
x=x0; y=y0; z=z0;
r_v_m(7,1:3)=[x,y,z];
x=x0; y=y0; z=z0+sz;
r_v_m(8,1:3)=[x,y,z];
for i=5:8
    r_v_m(i,4:5)=[1.5*pi,pi/2.0];
end

%face 3
x=x0; y=y0; z=z0+sz;
r_v_m(9,1:3)=[x,y,z];
x=x0; y=y0; z=z0;
r_v_m(10,1:3)=[x,y,z];
x=x0+sx; y=y0; z=z0;
r_v_m(11,1:3)=[x,y,z];
x=x0+sx; y=y0; z=z0+sz;
r_v_m(12,1:3)=[x,y,z];
for i=9:12
    r_v_m(i,4:5)=[pi,pi/2.0];
end

%face 4
x=x0+sx; y=y0; z=z0+sz;
r_v_m(13,1:3)=[x,y,z];
x=x0+sx; y=y0; z=z0;
r_v_m(14,1:3)=[x,y,z];
x=x0+sx; y=y0+sy; z=z0;
r_v_m(15,1:3)=[x,y,z];
x=x0+sx; y=y0+sy; z=z0+sz;
r_v_m(16,1:3)=[x,y,z];
for i=13:16
    r_v_m(i,4:5)=[pi/2.0,pi/2.0];
end

%face 5
x=x0+sx; y=y0+sy; z=z0+sz;
r_v_m(17,1:3)=[x,y,z];
x=x0+sx; y=y0+sy; z=z0;
r_v_m(18,1:3)=[x,y,z];
x=x0; y=y0+sy; z=z0;
r_v_m(19,1:3)=[x,y,z];
x=x0; y=y0+sy; z=z0+sz;
r_v_m(20,1:3)=[x,y,z];
for i=17:20
    r_v_m(i,4:5)=[0,pi/2.0];
end

%face 6, subsurface
x=x0+sx; y=y0; z=z0;
r_v_m(21,1:3)=[x,y,z];
x=x0; y=y0; z=z0;
r_v_m(22,1:3)=[x,y,z];
x=x0; y=y0+sy; z=z0;
r_v_m(23,1:3)=[x,y,z];
x=x0+sx; y=y0+sy; z=z0;
r_v_m(24,1:3)=[x,y,z];
for i=21:24
    r_v_m(i,4:5)=[pi/2.0,pi];
end

t_m=[]; %store coordinates of the points which is the intersected points between joint and cuboid in global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate the joint plane equation in global coordinate system%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=point_c;
alpha=pi/2-j_dd_dip(1,1);  %the rotato angle, trend of joint substract trend of x-axis
beta=j_dd_dip(2,1);  %the rotato angle of dip 
%transportation matrix, transport the local coordinate system into global coordinate system
nto=[cos(alpha)*cos(beta), -sin(alpha), cos(alpha)*sin(beta);
     sin(alpha)*cos(beta), cos(alpha),  sin(alpha)*sin(beta);
     -sin(beta),           0,           cos(beta);];
j_ov=nto*[0.0;0.0;10]; %normal vector of the joint plane in global coordinate system
o_lc=[j_ov(1,1);j_ov(2,1);j_ov(3,1);-(j_ov')*cs_o;];  %generate coefficients vector of the joint plane equation (oax+oby+ocz+od=0) in global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate the intersected point of infinite joint plane and cuboid face in loop%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:6;
	j=4*i-3; % to coincide with cuboid matrix
	cs_o=r_v_m(j,1:3)';  %attention, coordinates are stored in column
	alpha=pi/2-r_v_m(j,4);
	beta=r_v_m(j,5);
	nto=[cos(alpha)*cos(beta), -sin(alpha), cos(alpha)*sin(beta);
            sin(alpha)*cos(beta), cos(alpha),  sin(alpha)*sin(beta);
            -sin(beta),           0,           cos(beta);];
	%transformation matrix, transform the global coordinate system into local coordinate system
	otn=[cos(alpha)*cos(beta), sin(alpha)*cos(beta),  -sin(beta);
            -sin(alpha),          cos(alpha),            0;
            cos(alpha)*sin(beta), sin(alpha)*sin(beta),  cos(beta);];
	o_p_v=r_v_m(j:(j+3),1:3)';  %pop out the global coordinates of rectangular vertexes into o_p_v, pay attention to that the coordinates are stored in column
	n_p_v=otn*(o_p_v-[cs_o,cs_o,cs_o,cs_o]); %transform the global coordinates of rectangular vertexes into local coordinate 
	%generate coefficients vector of the joint plane equation
	%(nax+nby+ncz+nd=0) in local coordinate system
    n_lc=[((j_ov')*nto(:,1));((j_ov')*nto(:,2));((j_ov')*nto(:,3));((j_ov')*cs_o+o_lc(4,1));];
    %n_lc=[j_nv(1,1);j_nv(2,1);j_nv(3,1); -(otn*point_c)'*j_nv;]; 
	n_lc=PruneMartix(n_lc,1e-5);
	%interval of the cut window in local coordinate system
	xwmin=min(n_p_v(1,:)); xwmax=max(n_p_v(1,:));
	ywmin=min(n_p_v(2,:)); ywmax=max(n_p_v(2,:));
	%let z=0£¬then the euqation of the joint plane is a infinite line
	if (n_lc(1,1)==0)&&(n_lc(2,1)==0)  %joint plane parallel with the cuboid face
		continue
	end
	if n_lc(2,1)==0
		x1=-n_lc(4,1)/n_lc(1,1);  y1=ywmin;  x2=x1;  y2=ywmax;
	else
		x1=xwmin;  y1=(-n_lc(4,1)-n_lc(1,1)*x1)/n_lc(2,1);  
        x2=xwmax;  y2=(-n_lc(4,1)-n_lc(1,1)*x2)/n_lc(2,1);
	end
	%Liang youdong-Barsky intersection algorithm
	[p1,p2,is_flag]=LBLine2D(xwmin,ywmin,xwmax,ywmax,x1,y1,x2,y2);
	if is_flag     %if they are intersected, stored the globle coordinates of intersection point
   		np_j_v=[p1,p2];
		[row,col]=size(np_j_v);
		op_j_v=zeros(row,col);
			for k=1:col
				op_j_v(:,k)=nto*np_j_v(:,k)+cs_o;
			end
		 t_m=[t_m,op_j_v];     
	end     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%merge the same column of coordinates of the intersected points to obtain the vertexes of the joint%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_m=PruneMartix(t_m,1e-5);
[row,col]=size(t_m);
for m=1:col
	for n=(m+1):col
		if sum(abs((t_m(:,m)-t_m(:,n))))<1e-3
			t_m(:,n)=t_m(:,m);
		end
	end
end
o_j_v=unique((t_m'),'rows');
o_j_v=o_j_v';
[row,col]=size(o_j_v);
if col<=2   %at least 3 point, otherwise no aera of joint in the cuboid
	cv_flag=-1;
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sort the vertexes of the joint in joint local coordinate system to simulate joint in the concerned cuboid %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=point_c;
alpha=pi/2-j_dd_dip(1,1);
beta=j_dd_dip(2,1);
otn=[cos(alpha)*cos(beta), sin(alpha)*cos(beta), -sin(beta);
     -sin(alpha),          cos(alpha),           0;
     cos(alpha)*sin(beta), sin(alpha)*sin(beta), cos(beta);];
[row,col]=size(o_j_v);
t_m=zeros(row,col);
for m=1:col
	t_m(:,m)=o_j_v(:,m)-cs_o;
end
n_j_v=otn*t_m;
n_j_v=PruneMartix(n_j_v,1e-5);

o_j_p=get(h,'Vertices');
o_j_p=o_j_p';
[row,col]=size(o_j_p);
t_m=zeros(row,col);
for m=1:col
    t_m(:,m)=o_j_p(:,m)-cs_o;
end
n_j_p=otn*t_m;
n_j_p=PruneMartix(n_j_p,1e-5);


%%in joint local coordinate system, find the maximun interval of the joint
xwmin=min(n_j_v(1,:));
t_m=n_j_v(1,:);
ind=find(abs(n_j_v(1,:)-xwmin)<=1e-5);
t_m(ind)=xwmin;
n_j_v(1,:)=t_m;
t_m=n_j_v(2,:);
ywmin=min(t_m(ind)); 
xwmax=max(n_j_v(1,:)); 
t_m=n_j_v(1,:);
ind=find(abs(n_j_v(1,:)-xwmax)<=1e-5);
t_m(ind)=xwmax;
n_j_v(1,:)=t_m;
t_m=n_j_v(2,:);
ywmax=max(t_m(ind)); 

%determine the upper and lower part of vector AB
[row,col]=size(n_j_v);
point_u=[];   point_d=[];
for n=1:col  %sort the vertexes of the joint
	x1=n_j_v(1,n);   
	y1=n_j_v(2,n);
	ud_flag=(ywmin-ywmax)*x1+( xwmax-xwmin)*y1+(xwmin*ywmax-xwmax*ywmin);  %linear programming
	if abs(ud_flag)<1e-5
		ud_flag=0;
	end
	if ud_flag>0
		point_u=[point_u,n_j_v(:,n)];
	end
	if ud_flag<0
		point_d=[point_d,n_j_v(:,n)];
	end
end
[row,col]=size(point_u);
if col>1
	[~,ind]=sort(point_u(1,:),'ascend');
	point_u=point_u(:,ind);
end
[row,col]=size(point_d);
if col>1  
	[~,ind]=sort(point_d(1,:),'descend');
	point_d=point_d(:,ind);
end
n_j_v=[[xwmin;ywmin;0],point_u,[xwmax;ywmax;0],point_d,[xwmin;ywmin;0]];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%clear the useless varible%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear t_m;
clear point_u;
clear point_d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%n_j_v--the vertex coordinates of the polygon in local system which is formed by the infinite joint plane cutted by box 
%n_j_p--the vertex coordinates of the whole joint in local system

cv_flag=-1;
[row,col]=size(n_j_v);  pic_flag=0;
for i=1:col
	in_flag=IsPointInPoly2D(n_j_p,n_j_v(:,i));
	if (in_flag==1)||(in_flag==0)
		pic_flag=pic_flag+in_flag;
	end
	if pic_flag>=1
		cv_flag=1;
		return
	end	
end

[row,col]=size(n_j_p);  cip_flag=0;
for i=1:col
	in_flag=IsPointInPoly2D(n_j_v,n_j_p(:,i));
	if (in_flag==1)||(in_flag==0)
		cip_flag=cip_flag+in_flag;
	end
	if cip_flag>=1
		cv_flag=1;
		return
	end			
end

return

end