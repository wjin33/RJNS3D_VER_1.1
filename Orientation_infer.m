clear all; close all; clc;
t=load('采样裂隙倾向倾角.txt');  %j_dd_dip--2xn matrix, 节理的[倾向;倾角] 角度
t=t';
j_dd_dip=t*pi/180;
[jg_dd_dip,es_k]=EsFisherPara(j_dd_dip);
jg_dd_dip=jg_dd_dip*180/pi;
% %输出变量，每一个变量输出txt文件
jg_dd_dip %命名为“平均产状”
es_k      %命名为“Fisher 控制参数”



rmc=[10;10;10];%测窗的空间坐标，[x;y;z]
mwdy=5;%测窗半高
mwdx=5;%测窗半长
t=[135;90];%测窗产状   [倾向;倾角] 角度
mw_dd_dip=t*pi/180;
m=10;%Possion圆盘直径均值
w=CalWeightSamBias(j_dd_dip,rmc,mwdy,mwdx,mw_dd_dip,m);
[jg_dd_dip,es_k,flag]=EsFisherParaW(j_dd_dip,w);
jg_dd_dip=jg_dd_dip*180/pi;
% %输出变量，每一个变量输出txt文件
jg_dd_dip %命名为“平均产状”
es_k      %命名为“Fisher 控制参数”




es_lumtaa=0.12;%裂隙面积密度
es_lumtav=CalJointVolDens(jg_dd_dip,mw_dd_dip,es_lumtaa,m);
% %输出变量，每一个变量输出txt文件
es_lumtav%体积密度