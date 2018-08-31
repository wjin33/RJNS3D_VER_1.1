%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=0;  y0=0;  z0=0; %the converx of the Cuboid, i.e. the simulated area!
x1=80; y1=80; z1=80;
[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,x1,y1,z1);
DisAxis([0;0;0;],3); view(3);  axis equal


% Using SplineFit.m and the observed trace in Measure Windows to get the Spline PDF of trace which is approximated by 
% legendre orthogonal polynomial f(l), then using the obtained polynomial f(l) to calculate the PDF of Radius g(x);
function SplineFit.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oa=????
pfn=????
function [m,mc,g]=PoissonDiskDistr(oa,pfn,dmax,dstart)
%(pfn+1)--the oder of the function which need to be fitted
%dsmax--the maximum diameter 
%dstart--the initial diameter
%oa--the fitted coefficient the polynomial of trace, from low to high
%m,mc--A order moments, mean diameter of the joint disk
%g--diameter distribution function, it's a symbolic function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radius=DiskDistrRand(g,dsmax,dstart,num_joint)
%g--diameter distribution function, it's a symbolic variable
%dsmax--the maximum diameter
%dstart--the initial diameter
%radius--generated each Poisson disk radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_dd_dip=
function [jg_dd_dip,es_k,flag]=EsFisherPara(j_dd_dip)
%j_dd_dip--2xn matrix, trend and dip of the investigated joint
%jg_dd_dip--2x1, the mean trend and dip
%es_k--the estimated K value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function es_j_dd_dip=GenFisherRand(jg_dd_dip,es_k,num_joint)
%es_j_dd_dip--2xn, the generated theta and phi parameter  
%jg_dd_dip--2x1, the mean trend and dip of the investigated joint
%es_k--the estimated K value, pay attention to choose the appropriate k according to the number of sample
%num_joint--number the joints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mw_strike=
function jt_dd_dip=TransJointDDDips(mw_strike,es_j_dd_dip)
%mw_strike--1x1 measure window trend, unit: radian
%mw_dip--dip of measure window: 90 degree
%es_j_dd_dip--2xnthe generated theta and phi parameter
%jt_dd_dip--2xn the transformed trend and dip of joint group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mw_dd_dip=
es_lumtaa=???
function es_lumtav=CalJointVolDens(jg_dd_dip,mw_dd_dip,es_lumtaa,m)
%es_lumtav--volume density of the joint group 
%jg_dd_dip--2x1 mean strike of the joint group, unit: radian
%mw_dd_dip--2x1 mean strike of the MW, unit: radian
%es_lumtaa--area density of the joint group 
%m--mean diameter of the joint disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=x1-x0;
dy=y1-y0;
dz=z1-z0;
point_c=PoissonProcess3D2(es_lumtav,dx,dy,dz);
%dx--length  dy--width(north direction)  dz--high
%pay attention to the unit
%point_c--3xn matrix, the coordinate of center of the joint
[row,num_joint]=size(point_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patch_handle=zeros(1,num_joint);
for i=1:num_joint
    patch_handle(1,i)=GenDiskJoint(point_c(:,i),jt_dd_dip(:,i),radius(i),r_v_m);
end
%point_c--3x1 matrix, the center point of joint
%r_v_m--24x5 matrix, vertexes of cuboid and and attitude of its faces
%jt_dd_dip--2x1 matrix£¬trend and dip of the joint which is simulated
%j_v_m--coordinates of vertexes of the joint





postprocesser



function [hmw2,hmw3,theta,len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiRMWA2(x0,y0,z0,dx,dy,dz,joint_c,j_dd_dip,hm,rmwa,mwh,mwv,mw_dd_dip,ratio_h,ratio_v)
%(x0,y0,z0)--coordinate of the first vertex of cuboid
%(dx,dy,dz)--length,width,hight of the cuboid
%joint_c--3xn matrix, center of each joint
%j_dd_dip--2xn matrix, trend and dip of each joint
%hm--1xn matrix, handle of each joint
%rmwa--3xn centra of the MWs
%mwv--length of the vertical axis of inner measure window
%mwh--length of the horizontal axis of inner measure window
%mw_dd_dip--the strike of measure window
%ratio_h--the ratio of horizontal axis  between outer and inner MW
%ratio_v--the ratio of vertical axis between outer and inner MW

function [hmw2,hmw3,theta,len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiCMWA2(x0,y0,z0,dx,dy,dz,joint_c,j_dd_dip,hm,cmwa,cmra,mw_dd_dip,ratio_h,ratio_v)
%(x0,y0,z0)--coordinate of the first vertex of cuboid
%(dx,dy,dz)--length,width,hight of the cuboid
%joint_c--3xn matrix, center of each joint
%j_dd_dip--2xn matrix, trend and dip of each joint
%hm--1xn matrix, handle of each joint
%cmwa--3xn centra of the MWs
%cmra-1xn radius of the inner circle
%mw_dd_dip--the strike of measure window
%ratio_h--the ratio of horizontal axis over radius, 
%ratio_v--the ratio of vertical  axis over radius,

function [theta,len_c_tra,len_d_tra,len_t_tra,plane_p]=GenRMWA3(x0,y0,z0,dx,dy,dz,joint_c,j_dd_dip,hm,rmc,mwdy,mwdx,mw_dd_dip)
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

function rec_cv=CalFractalDim3D(sx0,sy0,sz0,mwh,mwv,ms,hm,cpm,occm,rm)
%rec_cv--1xn, the number of  box which cover the joints relat to the scale
%sx0,sy0,sz0--initial point of the model
%mwv--half length of the vertical axis of inner measure window
%mwh--half length of the horizontal axis of inner measure window
%ms--1xn measurement scale matrix 
%hm--1xn handle of each joint
%cpm--3xn center of each joint
%occm--2xn strike of each joint
%rm--1xn radius of each joint

function [rec_cv,rec_cc]=CalFractalDim2D(mwh,mwv,ms,ipcm)
%rec_cv--1xn, the number of  square which cover the joints relat to the scale
%rec_cc--1xn, the number of joints covered relate to scale
%mwv--half length of the vertical axis of  measure window
%mwh--half length of the horizontal axis of  measure window
%ms--1xn measurement scale matrix
%ipcm--6xn the coordinates of each joint endpoints











