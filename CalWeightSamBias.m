function w=CalWeightSamBias(j_dd_dip,rmc,mwdy,mwdx,mw_dd_dip,m)
%j_dd_dip--2xn matrix, trend and dip of joint group
%mw_dd_dip--2x1 matrix, trend and dip of MW
%rmc--3x1 matrix, center of measure window in global coordinate
%mwdy--half length of the vertical axis of measure window
%mwdx--half length of the horizontal axis of measure window
%m--the mean diameter of the Poisson disk
%w--1xn calculated weight matrix

%%

%transform to local MW coordinate
cs_o=rmc;
alpha=pi/2-mw_dd_dip(1,1);
beta=mw_dd_dip(2,1);
nto=[cos(alpha)*cos(beta),-sin(alpha),cos(alpha)*sin(beta);
     sin(alpha)*cos(beta),cos(alpha),sin(alpha)*sin(beta);
     -sin(beta),0,cos(beta);];
otn=[cos(alpha)*cos(beta),sin(alpha)*cos(beta),-sin(beta);
    -sin(alpha),cos(alpha),0;
    cos(alpha)*sin(beta),sin(alpha)*sin(beta),cos(beta);];
%angular point of measure windom in local coordinate system
np_cn1=[mwdx;-mwdy;0]; np_cn2=[-mwdx;mwdy;0];
%transform to global coordinate system
p_cn1=nto*np_cn1+cs_o; p_cn1=PruneMartix(p_cn1,1e-5);
p_cn2=nto*np_cn2+cs_o; p_cn2=PruneMartix(p_cn2,1e-5);
%the four angular point in local coordinate system
nrec_cn1=np_cn1; nrec_cn2=[np_cn1(1,1);np_cn2(2,1);0];
nrec_cn3=np_cn2; nrec_cn4=[np_cn2(1,1);np_cn1(2,1);0];
%transform the four angular point in local coordinate system 
%to global coordinate system, anticlockwise
orec_cn1=nto*nrec_cn1+cs_o; orec_cn2=nto*nrec_cn2+cs_o; 
orec_cn3=nto*nrec_cn3+cs_o; orec_cn4=nto*nrec_cn4+cs_o; 
orec_cn1=PruneMartix(orec_cn1,1e-5); orec_cn2=PruneMartix(orec_cn2,1e-5); 
orec_cn3=PruneMartix(orec_cn3,1e-5); orec_cn4=PruneMartix(orec_cn4,1e-5); 
%the vector of width and length of MW
r=orec_cn2-orec_cn1; s=orec_cn4-orec_cn1;
r=r/(r(1,1)^2+r(2,1)^2+r(3,1)^2)^0.5;
s=s/(s(1,1)^2+s(2,1)^2+s(3,1)^2)^0.5;
r=PruneMartix(r,1e-5); s=PruneMartix(s,1e-5);
%the normal vector of the MW: p2 
alpha=mw_dd_dip(1,1); beta=mw_dd_dip(2,1);
%%set the spherical coordinates parameter, theta and phi, refer to the
%%book by Prof.Yu Qingchun
theta=beta;
if alpha<=0.5*pi
	phi=0.5*pi-alpha;
else
	phi=2.5*pi-alpha;	
end
p2=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta);];
p2=PruneMartix(p2,1e-5);


%%

[row,col]=size(j_dd_dip);
w=zeros(1,col);%initialize weight matrix
for j=1:col
	%the normal vector of each joint: p1
	alpha=j_dd_dip(1,j); beta=j_dd_dip(2,j);
	%%set the spherical coordinates parameter, theta and phi, refer to the
    %%book by Prof.Yu Qingchun
	theta=beta;
	if alpha<=0.5*pi
		phi=0.5*pi-alpha;
	else
		phi=2.5*pi-alpha;	
	end
	p1=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta);]; p1=PruneMartix(p1,1e-5);
	%calculate ds
	ds=(2*mwdx*abs(dot(s,p1))+2*mwdy*abs(dot(r,p1)))/(1-(dot(p1,p2))^2)^0.5;
	%calculate vt
	vt=4*mwdx*mwdy*m*(1-(dot(p1,p2))^2)^0.5+0.25*pi*m^2*(2*mwdx*abs(dot(s,p1))+2*mwdy*abs(dot(r,p1)));
	w(j)=1.0/vt;
end
%%refer to the book by Prof.Yu Qingchun
k=col/sum(w);
w=k*w;
end