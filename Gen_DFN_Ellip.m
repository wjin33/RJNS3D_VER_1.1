clear all; close all;
x0=0;y0=0;z0=0;
dx=20;dy=40;dz=20;
%lumta--density
lumtav=0.005;
beta=pi/2;
k_axis=2;
jg_dd_dip=[135;45]*pi/180;
es_k=9;%Fisher

point_c=PoissonProcess3D2(lumtav,dx,dy,dz);
[row,num_joint]=size(point_c);

%%
%Fisher
j_dd_dip=GenFisherRand(jg_dd_dip,es_k,num_joint);
%j_dd_dip=[jg_dd_dip(1,1)*ones(1,num_joint);jg_dd_dip(2,1)*ones(1,num_joint)];

%%


B1=2;%
B2=12;%
major_axis=unifrnd(B1,B2,1,num_joint);   %长轴---均匀分布

% lamda=5;%指数分布--均值
% major_axis = exprnd(lamda,1,num_joint);  %长轴---指数分布
% 
% MU=15;%对数正太分布--均值
% SIGMA=1;%对数正太分布--标准差
% major_axis = normrnd(MU,SIGMA,1,num_joint); %长轴---正太分布
% 
% MU=3;%对数正太分布--均值
% SIGMA=0.1;%对数正太分布--标准差
% major_axis = lognrnd(MU,SIGMA,1,num_joint); %长轴---对数正太分布


% a_min=1.117;%长轴最小值
% a_max=10;%长轴最大值
% major_axis=JointSizeRand_TraUni(a_max,a_min,num_joint);  %长轴--迹长为均匀分布

% a_min=1;%长轴最小值
% a_max=50;%长轴最大值
% D=1.5;%迹长分维
% major_axis=JointSizeRand_TraFractal(D,a_min,a_max,num_joint);  %长轴--迹长为分形分布

% load Size_PDF
% major_axis=JointSizeRand(g,a_max,a_min,num_joint);  %长轴--迹长为多项式分布


%%

[r_v_m,hcuoid]=GenCuboid(x0,y0,z0,dx,dy,dz);
DisAxis([0;0;0;],3); % the global coordinate red--z; blue--y; green--x;
view(3);  axis equal
patch_handle=zeros(1,num_joint+1);
% the determined parameters of elliptical joints
for i=1:num_joint
patch_handle(1,i)=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,r_v_m);
end

break


%%%%%%%%%%增加单条裂隙%%%%%%%%%%%%%%%%%%%%
num_joint=num_joint+1;
j_dd_dip(:,num_joint)=[135;45]*pi/180;%裂隙产状
point_c(:,num_joint)=[10;20;10];%中心点坐标
major_axis(:,num_joint)=50;%长轴大小
patch_handle(1,i+1)=GenEllipiseJoint(point_c(:,num_joint),j_dd_dip(:,num_joint),major_axis(:,num_joint),k_axis,beta,r_v_m);
% 输出变量，每一个变量输出excel文件
% point_c %3 x n matrix，命名为“裂隙中心点坐标”
% j_dd_dip %2 x n matrix，命名为“裂隙产状”
% major_axis %1x n matrix，命名为“裂隙长轴尺寸”
%  xlswrite('C:\Users\Administrator\Desktop\程序\裂隙中心点坐标.xls',point_c);
%  xlswrite('C:\Users\Administrator\Desktop\程序\裂隙产状.xls',j_dd_dip);
%  xlswrite('C:\Users\Administrator\Desktop\程序\裂隙长轴尺寸.xls',major_axis);


% %
% %%%%%%%%%%%%%%矩形单测窗%%%%%%%%%%%%%%%%%%%%%
% rmc=[0.5*dx;0.5*dy;0.5*dz];% the center coordinates of measuring widow
% mw_dd_dip=[3.0*pi/2.0;1*pi/2.0]; % the orientation of measuring widow
% mwh=16;  mwv=8; % the vertical and horizontal length of measuring widow
% [beta_loc,theta,len_c_tra,len_d_tra,len_t_tra,plane_p,hmw3]=GenRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmc,mwh,mwv,mw_dd_dip);
% 测窗切割后输出变量，每一个变量输出txt文件
% plane_p %6 x n matrix，在测窗平面局部坐标系下的迹线两端的坐标，命名为“迹线端点坐标”
% len_c_tra %命名为“未删节迹长”
% len_d_tra %命名为“一端删节迹长”
% len_t_tra %命名为“两端删节迹长”


% %%%%%%%%%%%%%外方内圆双测窗%%%%%%%%%%%%%%%%%%%%
% cmwa=[0.5*dx;0.5*dy;0.5*dz]; %测窗中心点坐标
% mw_dd_dip=[pi/2.0;pi/2.0];   %测窗的倾向倾角
% ratio_h=1.5; %ratio_h--外测窗水平轴/内测窗半径之比
% ratio_v=2;   %ratio_v--外测窗垂直轴/内测窗半径之比
% cmra=5;      %测窗圆形半径
% [hmw2,hmw3,len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiCMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,cmwa,cmra,mw_dd_dip,ratio_h,ratio_v);
% %测窗切割后输出变量，每一个变量输出txt文件
% len_tc_tra %命名为“内测窗未删节迹长”
% len_td_tra %命名为“内测窗一端删节迹长-不包括测窗外长度”
% len_tcd_tra %命名为“内测窗一端删节迹长-包括测窗外长度”
% len_tct_tra %命名为“内测窗两端删节迹长”


% %%%%%%%%%%%%%%外方内方双测窗%%%%%%%%%%%%%%%%%%%%
% rmwa=[0.5*dx;0.5*dy;0.5*dz];  %测窗中心点坐标
% mw_dd_dip=[0*pi/2.0;pi/2.0];    %测窗的倾向倾角
% ratio_h=2.3; %ratio_h--外测窗/内测窗水平轴之比
% ratio_v=2;   %ratio_v--外测窗/内测窗垂直轴之比
% mwh=8;%内测窗水平轴
% mwv=6;%内测窗垂直轴
% [len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiRMWA2(x0,y0,z0,dx,dy,dz,point_c,j_dd_dip,patch_handle,rmwa,mwh,mwv,mw_dd_dip,ratio_h,ratio_v);
% % %测窗切割后输出变量，每一个变量输出txt文件
% % len_tc_tra %命名为“内测窗未删节迹长”
% % len_td_tra %命名为“内测窗一端删节迹长-不包括测窗外长度”
% % len_tcd_tra %命名为“内测窗一端删节迹长-包括测窗外长度”
% % len_tct_tra %命名为“内测窗两端删节迹长”

%%
% % %%%%%%%%%%%%%%%测窗迹线分维覆盖计算%%%%%%%%%%%%%%%%%%%%%
% %设定方格尺度向量
% ms=[4 2 1 0.5 0.25];%方格尺度数组
% [rec_cv,rec_cc]=CalFractalDim2D(mwh,mwv,ms,plane_p);
% %输出txt文件
% rec_cv %各个尺度下覆盖裂隙的方格数目，命名为“覆盖方格数”
% 
% %取芯操作，计算体积分维
% mwh=8;%取芯长/宽
% mwv=8;%取芯高度
% [r_v_m,hcuoid]=GenCuboid((0.5*(dx-2*mwh)),(0.5*(dy-2*mwv)),(0.5*(dy-2*mwv)),(2*mwh),(2*mwh),(2*mwh));
% view(3);  axis equal
% cpoint_c=[];cj_dd_dip=[];cradius=[];cpatch_handle=[];
% for i=1:num_joint
%     h=GenEllipiseJoint(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,r_v_m);%记录与取芯区域相交的圆盘个数
%     if h~=-1
%         cpatch_handle=[cpatch_handle,h];
%         cpoint_c=[cpoint_c,point_c(:,i)];
%         cj_dd_dip=[cj_dd_dip,j_dd_dip(:,i)];
%         cradius=[cradius,major_axis(i)];
%     end
% end
% ms=[2.5 2 1.5 1 0.5];%方格尺度数组
% sx0=(0.5*(dx-2*mwh));  sy0=(0.5*(dy-2*mwv)); sz0=(0.5*(dy-2*mwv));
% rec_cv=CalFractalDim3D(sx0,sy0,sz0,mwh,mwv,ms,cpatch_handle,cpoint_c,cj_dd_dip,cradius);
% %输出txt文件
% rec_cv%各个尺度下覆盖裂隙的方格数目，命名为“三维覆盖方格数”
% 



%%%%%%%%%%%%%%%生成多个平行钻孔%%%%%%%%%%%%%%%%%%%%%
% num_borehole=5;%钻孔的个数，固定值
% m=20;%m-圆周方向的等分数，固定值
% 
% location=[10,5,10,15,10;
%           10,15,20,25,30;
%           0,0,0,0,0;];    % 钻孔底部空间坐标，需要全部输入
%       
% bore_dip_direction=1*pi/4;  % 钻孔倾向，需要输入
% bore_dip_angle=0;           % 钻孔倾角，需要输入
% bore_dip=[bore_dip_direction*ones(1,num_borehole);bore_dip_angle*ones(1,num_borehole);];
% H=[20 20 20 20 20];%各个钻孔的深度，需要全部输入
% radius=[1 1 1 1 1];%各个钻孔的半径，需要全部输入
% for j=1:num_borehole
%     Genborehole(location(:,j),bore_dip(:,j),H(:,j),radius(:,j),m);
% end
% 
% flag_int=zeros(num_joint,num_borehole);
% flag_full=zeros(num_joint,num_borehole);
% for i=1:num_joint
%     for j=1:num_borehole
%         [flag_int(i,j),flag_full(i,j)]=IsborejointInsect(point_c(:,i),j_dd_dip(:,i),major_axis(i),k_axis,beta,location(:,j),bore_dip(:,j),radius(:,j));
%     end
% end
% N_intersect=sum(flag_int);
% N_fullintsect=sum(flag_full);
% mean_size=zeros(num_borehole,1);
% sigma=zeros(num_borehole,1);
% for i=1:num_borehole
%     [mean_size(i,1),sigma(i,1)]=CalmeanvarSize(lumtav,j_dd_dip(:,1),k_axis,beta,bore_dip(:,i),radius(:,i),H(:,i),N_intersect(1,i),N_fullintsect(1,i));
% end
% sigma=mean(sigma);
% mean_size=mean(mean_size);
%%输出txt文件
%N_intersect %命名为“裂隙与各钻孔交切条数”
%N_fullintsect %命名为“裂隙与各钻孔完全交切条数”
%mean_size %命名为“推断裂隙尺寸均值”
%sigma  %命名为“推断裂隙尺寸标准差”


