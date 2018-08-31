function [np1,np2,is_jmw_flag]=IsElliJointMWInsect(joint_c,j_dd_dip,h_j,p_cn1,p_cn2,mw_dd_dip)
%np1,np2--3x1 matrix, the local coordinates of the endpoints of trace in Measuring window
%op1,op2--3x1 matrix, the global coordinates of the endpoints of trace in Measuring window
%p_cn1--3x1 matrix, the global coordinates of left lower point of measuring window
%p_cn2--3x1 matrix, the global coordinates of  right upper point of measuring window 
%np_cn1--3x1 matrix, the local coordinates of left lower point of measuring window
%np_cn2--3x1 matrix, the local coordinates of right upper point of measuring window 
%mw_dd_dip--2x1 matrix, the trend and dip of measuring window 
%j_dd_dip--2x1 matrix, the trend and dip of joint
%joint_c--3x1 matrix, the center point of joint
%h_j--handle of joint
%is_jmw_flag--1:intersect, 0:non-intersect


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate the measuring window equation in global coordinate system%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=p_cn1;
alpha=pi/2-mw_dd_dip(1,1);%the rotato angle, trend of measuring window substract trend of x-axis
phi=mw_dd_dip(2,1);
nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
     sin(alpha)*cos(phi),  cos(alpha), sin(alpha)*sin(phi);
        -sin(phi),             0,         cos(phi);];
wm_nv=nto*[0.0;0.0;10]; %normal vector of the measuring window in global coordinate system
o_lc=[wm_nv(1,1);wm_nv(2,1);wm_nv(3,1);-(wm_nv')*cs_o;];%generate coefficients vector of the measuring window equation (oax+oby+ocz+od=0) in global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=joint_c;
alpha=pi/2-j_dd_dip(1,1);
phi=j_dd_dip(2,1);
nto=[cos(alpha)*cos(phi), -sin(alpha), cos(alpha)*sin(phi);
     sin(alpha)*cos(phi), cos(alpha),  sin(alpha)*sin(phi);
     -sin(phi),           0,           cos(phi);];
otn=[ cos(alpha)*cos(phi), sin(alpha)*cos(phi), -sin(phi);
      -sin(alpha),          cos(alpha),          0;
      cos(alpha)*sin(phi), sin(alpha)*sin(phi), cos(phi);];
%generate coefficients vector of the measuring window equation (nax+nby+ncz+nd=0) in local joint coordinate system
n_lc=[((wm_nv')*nto(:,1));((wm_nv')*nto(:,2));((wm_nv')*nto(:,3));(o_lc(4,1)+(wm_nv')*cs_o);];
n_lc=PruneMartix(n_lc,1e-5);
ovm=get(h_j,'Vertices');
ovm=ovm';
[row,col]=size(ovm);
nvm=zeros(row,col);
for i=1:col
	nvm(:,i)=otn*(ovm(:,i)-cs_o);
end
nvm=PruneMartix(nvm,1e-5);
%%in joint local coordinate system, find the maximun interval of the joint
xwmin=min(nvm(1,:)); xwmax=max(nvm(1,:)); ywmin=min(nvm(2,:)); ywmax=max(nvm(2,:)); 
%let z=0, then the euqation of the measuring window  plane is a infinite line
if (n_lc(2,1)==0)&&(n_lc(1,1)==0)  %joint plane parallel with  measuring window 
	is_jmw_flag=0;
	return
end
if n_lc(2,1)==0
	x1=-n_lc(4,1)/n_lc(1,1);  y1=ywmin;  x2=x1;  y2=ywmax;
else
	x1=xwmin; y1=(-n_lc(4,1)-n_lc(1,1)*x1)/n_lc(2,1);  x2=xwmax;  y2=(-n_lc(4,1)-n_lc(1,1)*x2)/n_lc(2,1);
end
npl=[]; %store the point of intersection between polygon and line segment
%Liang youdong-Barsky intersection algorithm
[np1,np2,is_flag]=LBLine2D(xwmin,ywmin,xwmax,ywmax,x1,y1,x2,y2);   %if they are intersected, stored the local coordinates of intersection point
if is_flag==0
	is_jmw_flag=0;
	return
else
	nvm=[nvm, nvm(:,1)];
	[row,col]=size(nvm);
	for i=1:(col-1)
		poly_p1=nvm(:,i);
		poly_p2=nvm(:,(i+1));
		[insc_p,is_ls_flag]=GetLSIP2D(np1,np2,poly_p1,poly_p2);
		if is_ls_flag~=0
			npl=[npl,insc_p];
		end
	end
	[row,col]=size(npl);
	if col~=2   %tangent to the polygon
		is_jmw_flag=0;
		return
	end
end  
opl=zeros(row,col);
for i=1:col
	opl(:,i)=nto*npl(:,i)+cs_o;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs_o=p_cn1;
alpha=pi/2-mw_dd_dip(1,1);
phi=mw_dd_dip(2,1);
nto=[cos(alpha)*cos(phi),-sin(alpha),cos(alpha)*sin(phi);sin(alpha)*cos(phi),cos(alpha),sin(alpha)*sin(phi);-sin(phi),0,cos(phi);];
otn=[cos(alpha)*cos(phi),sin(alpha)*cos(phi),-sin(phi);-sin(alpha),cos(alpha),0;cos(alpha)*sin(phi),sin(alpha)*sin(phi),cos(phi);];
np_cn1=otn*(p_cn1-cs_o);  np_cn1=PruneMartix(np_cn1,1e-5); 
np_cn2=otn*(p_cn2-cs_o);  np_cn2=PruneMartix(np_cn2,1e-5);
npl(:,1)=otn*(opl(:,1)-cs_o); 
npl(:,2)=otn*(opl(:,2)-cs_o);
xwmin=min([np_cn1(1,1),np_cn2(1,1)]); xwmax=max([np_cn1(1,1),np_cn2(1,1)]);
ywmin=min([np_cn1(2,1),np_cn2(2,1)]); ywmax=max([np_cn1(2,1),np_cn2(2,1)]);
x1=npl(1,1); y1=npl(2,1); x2=npl(1,2); y2=npl(2,2);
[np1,np2,is_flag]=LBLine2D(xwmin,ywmin,xwmax,ywmax,x1,y1,x2,y2);
if is_flag==0
	is_jmw_flag=0;
	return
else
	is_jmw_flag=1;
	np1=PruneMartix(np1,1e-5);  np2=PruneMartix(np2,1e-5); 
	op1=nto*np1+cs_o;   op2=nto*np2+cs_o;
	op1=PruneMartix(op1,1e-5);  op2=PruneMartix(op2,1e-5); 
	line(([op1(1,1),op2(1,1)]),([op1(2,1),op2(2,1)]),([op1(3,1),op2(3,1)]));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%