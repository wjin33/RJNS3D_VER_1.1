clear all; close all;
x0=0;y0=0;z0=0;
dx=20;dy=40;dz=20;
lumtav=0.005;
beta=pi/180*90;
k_axis=2;
jg_dd_dip=[135;45]*pi/180;
es_k=9;
point_c=PoissonProcess3D2(lumtav,dx,dy,dz);
[row,num_joint]=size(point_c);
% Fisher分布
 j_dd_dip=GenFisherRand(jg_dd_dip,es_k,num_joint);

%平行分布
%A2 j_dd_dip=[jg_dd_dip(1,1)*ones(1,num_joint);jg_dd_dip(2,1)*ones(1,num_joint)];

%B1 B1=@B1;
%B1 B2=@B2;
%B1 major_axis=unifrnd(B1,B2,num_joint);   %长轴---均匀分布

 lamda=5;
 major_axis = exprnd(lamda,1,num_joint);  %长轴---指数分布
% 
%B3 MU=@MU;
%B3 SIGMA=@SIGMA;
%B3 major_axis = normrnd(MU,SIGMA,1,num_joint); %长轴---正太分布

%B4 MU=@MU;
%B4 SIGMA=@SIGMA;
%B4 major_axis = lognrnd(MU,SIGMA,1,num_joint); %长轴---对数正太分布

%B5 a_min=@a_min;
%B5 a_max=@a_max;
%B5 major_axis=JointSizeRand_TraUni(a_max,a_min,num_joint);  %长轴--迹长为均匀分布

%B6 a_min=@a_min;
%B6 a_max=@a_max;
%B6 a_D=@a_D;
%B6 major_axis=JointSizeRand_TraFractal(a_D,a_min,a_max,num_joint);  %长轴--迹长为分形分布

[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,dx,dy,dz);
DisAxis([0;0;0;],3); 
view(3);  axis equal
patch_handle=zeros(1,num_joint);
for i=1:num_joint
patch_handle(1,i)=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,r_v_m);
end

%%%%%%%%%%增加单条裂隙%%%%%%%%%%%%%%%%%%%%
num_joint=num_joint+1;
j_dd_dip(:,num_joint)=[135;45]*pi/180;%裂隙产状
point_c(:,num_joint)=[10;10;10];%中心点坐标
major_axis(:,num_joint)=30;%长轴大小
patch_handle(1,i+1)=GenEllipiseJoint(point_c(:,num_joint),j_dd_dip(:,num_joint),major_axis(:,num_joint),k_axis,beta,r_v_m);
% 输出变量，每一个变量输出excel文件
% point_c %3 x n matrix，命名为“裂隙中心点坐标”
% j_dd_dip %2 x n matrix，命名为“裂隙产状”
% major_axis %1x n matrix，命名为“裂隙长轴尺寸”




%矩形单测窗
%C1 rmc=[0.5*dx;0.5*dy;0.5*dz];
%C1 mw_dd_dip=[@C1_T*pi/180;@C1_D*pi/180];
%C1 mwh=@C1_mwh;  mwv=@C1_mwv;
%C1 [beta_loc,theta,len_c_tra,len_d_tra,len_t_tra,plane_p,hmw3]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);
%C1 dlmwrite('@C1_path1',plane_p);
%C1 dlmwrite('@C1_path2',len_c_tra);
%C1 dlmwrite('@C1_path3',len_d_tra);
%C1 dlmwrite('@C1_path4',len_t_tra);

%外方内圆双测窗
%C2 cmwa=[0.5*dx;0.5*dy;0.5*dz];
%C2 mw_dd_dip=[@C2_T/180*pi;pi*@C2_D/180];
%C2 ratio_h=@C2_ratio_h; 
%C2 ratio_v=@C2_ratio_v;
%C2 cmra=@C2_cmra;
%C2 [hmw2,hmw3,len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiCMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,cmwa,cmra,mw_dd_dip,ratio_h,ratio_v);
%C2 dlmwrite('@C2_path1',len_tc_tra);
%C2 dlmwrite('@C2_path2',len_td_tra);
%C2 dlmwrite('@C2_path3',len_tcd_tra);
%C2 dlmwrite('@C2_path4',len_tct_tra);

%外方内方双测窗
%C3 rmwa=[0.5*dx;0.5*dy;0.5*dz];
%C3 mw_dd_dip=[@C3_T*pi/180;pi/180*@C3_D];
%C3 ratio_h=@C3_ratio_h;
%C3 ratio_v=@C3_ratio_v;
%C3 mwh=@C3_mwh;
%C3 mwv=@C3_mwv;
%C3 [len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmwa,mwh,mwv,mw_dd_dip,ratio_h,ratio_v);
%C3 dlmwrite('@C3_path1',len_tc_tra);
%C3 dlmwrite('@C3_path2',len_td_tra);
%C3 dlmwrite('@C3_path3',len_tcd_tra);
%C3 dlmwrite('@C3_path4',len_tct_tra);

%测窗迹线分维覆盖计算
%D ms=[@D_ms];
%D [rec_cv,rec_cc]=CalFractalDim2D(mwh,mwv,ms,plane_p);
%D dlmwrite('@D_path',rec_cv);

%计算体积分维
 mwh=8;
 mwv=8;
 [r_v_m,hcuoid]=GenCuboid((0.5*(dx-2*mwh)),(0.5*(dy-2*mwv)),(0.5*(dy-2*mwv)),(2*mwh),(2*mwh),(2*mwh));
 view(3);axis equal;
 cpoint_c=[];cj_dd_dip=[];cradius=[];cpatch_handle=[];
 for i=1:num_joint
     h=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,r_v_m);
     if h~=-1
         cpatch_handle=[cpatch_handle,h];
         cpoint_c=[cpoint_c,point_c(:,i)];
         cj_dd_dip=[cj_dd_dip,j_dd_dip(:,i)];
         cradius=[cradius,major_axis(i)];
    end
 end
 ms=[2.5,2,1.5,1,0.5];
 sx0=(0.5*(dx-2*mwh));  sy0=(0.5*(dy-2*mwv)); sz0=(0.5*(dy-2*mwv));
 rec_cv=CalFractalDim3D(sx0,sy0,sz0,mwh,mwv,ms,cpatch_handle,cpoint_c,cj_dd_dip,cradius);
 dlmwrite('D:\Graduate Work\封装工具箱\程序\椭圆\三维覆盖方格数.txt',rec_cv);

%生成多个平行钻孔
%F num_borehole=5;
%F m=@Fm;
%F location=[@Fx;
%F           @Fy;
%F           @Fz;];
%F bore_dip_direction=@F2*pi/180;
%F bore_dip_angle=@F1*pi/180;
%F bore_dip=[bore_dip_direction*ones(1,num_borehole);bore_dip_angle*ones(1,num_borehole);];
%F H=[@FH];
%F radius=[@FR];
%F for j=1:num_borehole
%F     Genborehole(location(:,j),bore_dip(:,j),H(:,j),radius(:,j),m);
%F end
%F flag_int=zeros(num_joint,num_borehole);
%F flag_full=zeros(num_joint,num_borehole);
%F for i=1:num_joint
%F     for j=1:num_borehole
%F         [flag_int(i,j),flag_full(i,j)]=IsborejointInsect(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,location(:,j),bore_dip(:,j),radius(:,j));
%F     end
%F end
%F N_intersect=sum(flag_int);
%F N_fullintsect=sum(flag_full);
%F mean_size=zeros(num_borehole,1);
%F sigma=zeros(num_borehole,1);
%F for i=1:num_borehole
%F     [mean_size(i,1),sigma(i,1)]=CalmeanvarSize(lumtav,j_dd_dip(:,1),k_axis,beta,bore_dip(:,i),radius(:,i),H(:,i),N_intersect(1,i),N_fullintsect(1,i));
%F end
%F sigma=mean(sigma);
%F mean_size=mean(mean_size);
%F dlmwrite('@F_path1',N_intersect);
%F dlmwrite('@F_path2',N_fullintsect);
%F dlmwrite('@F_path3',mean_size);
%F dlmwrite('@F_path4',sigma);
helpdlg('计算完毕');


