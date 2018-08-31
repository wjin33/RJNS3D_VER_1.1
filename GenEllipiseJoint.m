function h=GenEllipiseJoint(point_c,j_dd_dip,major_axis,k_axis,beta,r_v_m)
%h--handle of the figure£¬when h=-1£¬means non-intersect
%point_c--3x1 matrix, the center point of joint
%r_v_m--24x5 matrix, vertexes of cuboid and and attitude of its faces
%cs_o--3x1 matrix, coordinate of the original point
%j_dd_dip--2x1 matrix£¬trend and dip of the joint which is simulated
%n_j_v--[] matrix, coordinates of vertexes of the joint in local coordinate system
%o_j_v--[] matrix, coordinates of vertexes of the joint in global coordinate system
%np_j_v--[] matrix, coordinates of vertexes of the joint in local coordinate system of the cut window
%op_j_v--[] matrix, coordinates of vertexes of the joint in global coordinate system of the cut window
%n_lc--coefficients vector of the joint plane equation in local coordinate system
%o_lc--coefficients vector of the joint plane equation in global coordinate system
%t_m--[] matrix, for temporarily using
%j_v_m--coordinates of vertexes of the joint
%major_axis--the characteristic dimension
%k_axis--the ratio of major to minor axis of a ellipse
%beta--the rotation angle
%cp_nvm--subset of vertexes between polygon and circle

t_m=[]; %store coordinates of the points which is the intersected points between joint and cuboid in global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate the joint plane equation in global coordinate system%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=point_c;
alpha=pi/2-j_dd_dip(1,1);  %the rotato angle, trend of joint substract trend of x-axis
phi=j_dd_dip(2,1);  %the rotato angle of dip 
%transportation matrix, transport the local coordinate system into global coordinate system
nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
     sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
     -sin(phi),           0,           cos(phi);];
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
	phi=r_v_m(j,5);
	nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
            sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
            -sin(phi),           0,           cos(phi);];
	%transformation matrix, transform the global coordinate system into local coordinate system
	otn=[cos(alpha)*cos(phi), sin(alpha)*cos(phi),  -sin(phi);
            -sin(alpha),          cos(alpha),            0;
            cos(alpha)*sin(phi), sin(alpha)*sin(phi),  cos(phi);];
	o_p_v=r_v_m(j:(j+3),1:3)';  %pop out the global coordinates of rectangular vertexes into o_p_v, pay attention to that the coordinates are stored in column
	n_p_v=otn*(o_p_v-[cs_o,cs_o,cs_o,cs_o]); %transform the global coordinates of rectangular vertexes into local coordinate 
	%generate coefficients vector of the joint plane equation
	%(nax+nby+ncz+nd=0) in local coordinate system
	j_nv=otn*j_ov;
	n_lc=[j_nv(1,1);j_nv(2,1);j_nv(3,1); -(otn*(point_c-cs_o))'*j_nv;]; 
	%n_lc=[((j_ov')*nto(:,1));((j_ov')*nto(:,2));((j_ov')*nto(:,3));((j_ov')*cs_o+o_lc(4,1));];
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
		x1=xwmin; y1=(-n_lc(4,1)-n_lc(1,1)*x1)/n_lc(2,1);  
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
	h=-1;
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%clear the useless varible%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear t_m;
clear op_j_v;
clear np_j_v
clear o_p_v;
clear n_p_v;
clear r_v_m;
clear n_lc;
clear o_lc;
clear j_ov;
clear j_nv;
% pack   PACK performs memory garbage collection. Extended MATLAB
%     sessions may cause memory to become fragmented, preventing
%     large variables from being stored. PACK is a command that
%     saves all variables on disk, clears the memory, and then
%     reloads the variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sort the vertexes of the joint in joint local coordinate system to simulate joint in the concerned cuboid %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=point_c;
alpha=pi/2-j_dd_dip(1,1);
phi=j_dd_dip(2,1);
otn=[cos(alpha)*cos(phi), sin(alpha)*cos(phi), -sin(phi);
     -sin(alpha),          cos(alpha),           0;
     cos(alpha)*sin(phi), sin(alpha)*sin(phi), cos(phi);];

nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
     sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
     -sin(phi),           0,           cos(phi);];
[row,col]=size(o_j_v);
t_m=zeros(row,col);
for m=1:col
	t_m(:,m)=o_j_v(:,m)-cs_o;
end
n_j_v=otn*t_m;
n_j_v=PruneMartix(n_j_v,1e-5);

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
	[t_m,ind]=sort(point_u(1,:),'ascend');
	point_u=point_u(:,ind);
end
[row,col]=size(point_d);
if col>1  
	[t_m,ind]=sort(point_d(1,:),'descend');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       condition 1: ellipise is inside the polygon                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%count the number of vertexes of circle which is in the joint polygon
eip_nvm=[];
%calculate the point of the ellipise trace
    Beta = beta;
    sinbeta = sin(Beta);
    cosbeta = cos(Beta);
    angle = linspace(0, 360, 30).* (pi / 180);
    sinalpha = sin(angle);
    cosalpha = cos(angle);
    [row,col]=size(angle);
    e_nvm=zeros(3,col);
    for i=1:col
        e_nvm(1,i) = major_axis/2 * cosalpha(1,i) * cosbeta - major_axis/(2*k_axis) * sinalpha(1,i) * sinbeta;
        e_nvm(2,i) = major_axis/2 * cosalpha(1,i) * sinbeta + major_axis/(2*k_axis) * sinalpha(1,i) * cosbeta;
    end
    e_nvm=PruneMartix(e_nvm,1e-5);
    
% figure(2)
% colort=rand(1,3);
% h=patch('Vertices',(e_nvm'),'Faces',(1:1:col),'FaceColor',colort,'EdgeColor',colort);
% hold on
% plot(n_j_v(1,:),n_j_v(2,:),'k-o','LineWidth',2,'MarkerSize',10)


      
%determine the vertexes of ellipise inside the polygon or not
for i=1:col 
	in_flag=IsPointInPoly2D(n_j_v,e_nvm(:,i));
	if (in_flag==1)||(in_flag==0)
		eip_nvm=[eip_nvm,e_nvm(:,i)];
	end
end
[row_eip,col_eip]=size(eip_nvm);
if col_eip==col
	n_j_v=eip_nvm;
	[row,col]=size(n_j_v);
	t_m=nto*n_j_v;
	for m=1:col
		o_j_v(:,m)=t_m(:,m)+cs_o;
	end
	o_j_v=PruneMartix(o_j_v,1e-5);
	colort=rand(1,3);
	h=patch('Vertices',(o_j_v'),'Faces',(1:1:col),'FaceColor',colort,'EdgeColor',colort);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%clear the useless varible%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	clear e_nvm;
	clear t_m;
	clear n_j_v;
	clear o_j_v;
	clear eip_nvm;
    clear row;
    clear row_eip;
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            condition 2: polygon is inside the ellipise           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine the vertexes of ellipise inside the polygon or not
pie_nvm=[];
[row,col]=size(n_j_v);
for i=1:col 
	in_flag=IsPointInPoly2D(e_nvm,n_j_v(:,i));
	if (in_flag==1)||(in_flag==0)
		pie_nvm=[pie_nvm,n_j_v(:,i)];
	end
end
[row_pie,col_pie]=size(pie_nvm);
if col_pie==col
	n_j_v=pie_nvm;
	[row,col]=size(n_j_v);
	t_m=nto*n_j_v;
	for m=1:col
		o_j_v(:,m)=t_m(:,m)+cs_o;
	end
	o_j_v=PruneMartix(o_j_v,1e-5);
	colort=rand(1,3);
	h=patch('Vertices',(o_j_v'),'Faces',(1:1:col),'FaceColor',colort,'EdgeColor',colort);
	%%%clear the useless varible
	clear t_m;
	%clear n_j_v;
	clear o_j_v;
	clear e_nvm;
	clear pie_nvm;
    clear row_pie;
    clear row;
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 condition 3: intersection                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epi_nvm=[];
[row,col]=size(n_j_v); 
	OTN=[cos(beta), -sin(beta),0; 
         sin(beta),  cos(beta),0;
             0    ,      0    ,1;];
	%transformation matrix, transform the global coordinate system into local coordinate system
	NTO=[ cos(beta),  sin(beta),0; 
         -sin(beta),  cos(beta),0;
              0    ,      0    ,1;];
% plot(OTN*n_j_v (1,:),OTN*n_j_v (2,:),'r--*','LineWidth',2,'MarkerSize',10)
         
for i=1:(col-1)  %%find the edges of polygon which intersect with ellipise 
	poly_p1=NTO*n_j_v(:,i);
	poly_p2=NTO*n_j_v(:,(i+1));
	[insc_p,is_flag]=GetELSIP2D([0;0;0],major_axis,k_axis,poly_p1,poly_p2);
	if is_flag~=0
		epi_nvm=[epi_nvm,OTN*insc_p];
	end
end

% hold on
% plot(eip_nvm(1,:),eip_nvm(2,:),'r--*','LineWidth',2,'MarkerSize',10)
% plot(pie_nvm(1,:),pie_nvm(2,:),'r--*','LineWidth',2,'MarkerSize',10)
%  figure(3)
%  m=NTO*e_nvm;
%  plot(m(1,:),m(2,:),'r--','LineWidth',2,'MarkerSize',10)
%  n=NTO*n_j_v;
%  plot(n(1,:),n(2,:),'r-','LineWidth',2,'MarkerSize',10)
%  rectangle('position',[-major_axis/2,-major_axis/(2*k_axis),major_axis,major_axis/(k_axis)],'curvature',[1,1])
%  hold on
% plot(epi_nvm(1,:),epi_nvm(2,:),'rO','LineWidth',2,'MarkerSize',10)
% epi_nvm=[epi_nvm(1,:);
%          epi_nvm(2,:);
%          zeros(size(epi_nvm(2,:)));];


cp_nvm=[eip_nvm,pie_nvm,epi_nvm];
cp_nvm=PruneMartix(cp_nvm,1e-5);
%%%clear the useless varible
clear e_nvm;
clear n_j_v;
clear o_j_v;
clear poly_p1;
clear poly_p2;
clear insc_p;
clear eip_nvm;
clear pie_nvm;
clear epi_nvm;
[~,col]=size(cp_nvm);
for m=1:col
	for n=(m+1):col
		if sum(abs((cp_nvm(:,m)-cp_nvm(:,n))))<1e-3 
		cp_nvm(:,n)=cp_nvm(:,m);
		end
	end
end
cp_nvm=PruneMartix(cp_nvm,1e-5);
cp_nvm=unique((cp_nvm'),'rows');
cp_nvm=cp_nvm';
[~,col]=size(cp_nvm);
if col<=2
	h=-1;
	return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sort the point of cp_nvm in clockwise
xwmin=min(cp_nvm(1,:));
t_m=cp_nvm(1,:);
ind=find(abs(cp_nvm(1,:)-xwmin)<=1e-5);
t_m(ind)=xwmin;
cp_nvm(1,:)=t_m;
t_m=cp_nvm(2,:);
ywmin=min(t_m(ind)); 
xwmax=max(cp_nvm(1,:)); 
t_m=cp_nvm(1,:);
ind=find(abs(cp_nvm(1,:)-xwmax)<=1e-5);
t_m(ind)=xwmax;
cp_nvm(1,:)=t_m;
t_m=cp_nvm(2,:);
ywmax=max(t_m(ind)); 
[~,col]=size(cp_nvm);
point_u=[];   
point_d=[];
for n=1:col
	x1=cp_nvm(1,n);   y1=cp_nvm(2,n);
	ud_flag=(ywmin-ywmax)*x1+( xwmax-xwmin)*y1+(xwmin*ywmax-xwmax*ywmin);
	if abs(ud_flag)<1e-5
		ud_flag=0;
	end
	if ud_flag>0
		point_u=[point_u,cp_nvm(:,n)];
	end
	if ud_flag<0
		point_d=[point_d,cp_nvm(:,n)];
	end
end
[row,col]=size(point_u);
if col>1
	[t_m,ind]=sort(point_u(1,:),'ascend');
	point_u=point_u(:,ind);
end
[~,col]=size(point_d);
if col>1
	[~,ind]=sort(point_d(1,:),'descend');
	point_d=point_d(:,ind);
end
cp_nvm=[[xwmin;ywmin;0],point_u,[xwmax;ywmax;0],point_d,[xwmin;ywmin;0]];
n_j_v=cp_nvm;
[~,col]=size(n_j_v);
t_m=nto*n_j_v;
for m=1:col
	o_j_v(:,m)=t_m(:,m)+cs_o;
end
o_j_v=PruneMartix(o_j_v,1e-5);
colort=rand(1,3);
h=patch('Vertices',(o_j_v'),'Faces',(1:1:col),'FaceColor',colort,'EdgeColor',colort);
%%%clear the useless varible
clear t_m;
clear n_j_v;
clear o_j_v;
clear point_u;
clear point_d;
clear cp_nvm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

