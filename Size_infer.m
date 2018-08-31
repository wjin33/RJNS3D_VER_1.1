close all; clear all; clc
N_c_tra=1321;%未删节迹线条数
N_d_tra=508;%一端未删节迹线条数
mwh=25;%测窗水平半长
mwv=18;%测窗竖直半高
theta=0.8861;%迹线与测窗水平线夹角，弧度
len_c_tra=load('未删节迹长数据.txt'); %len_c_tra--1xn matrix, 未删节迹长列表
Numb_hist=20;%直方个数
opt_n=8;  %迹长拟合阶数


dens_aera= (2*N_c_tra+N_d_tra)/(2*4*mwh*mwv);
L=[];PDF=[];
for i=3:Numb_hist
    [F,X] = ecdf(len_c_tra,'Function','cdf');  % compute empirical cumulative distribution function cdf
    Bin_.rule = 3;
    Bin_.nbins = i;
    [C,E] = dfswitchyard('dfhistbins',len_c_tra,[],[],Bin_,F,X);
    [N,C] = ecdfhist(F,X,'edges',E); % empirical pdf from cdf
    nc=N*length(len_c_tra);
    dh=((nc./dens_aera)./(cos(theta)*sin(theta).*(C.^2)-(2*mwh*sin(theta)+2*mwv*cos(theta)).*C+4*mwh*mwv));
    L=[L,C]; PDF=[PDF,dh];
end
figure
bar(C,nc/sum(N),'w')
box off
xlabel('trace length')
ylabel('统计直方')
[bm1,ind] = sort(L(1,:),'ascend');
L=L(1,ind);PDF=PDF(1,ind);
l_max=max(len_c_tra);
l_min=min(len_c_tra);

[oa,res,aic]=lsqfun2(opt_n,L,PDF,l_max);
	str_f=strcat('x.*(x-',num2str(l_max),').*(');
	for i=1:(opt_n-1)
	str_f=strcat(str_f,num2str(oa(i),'%+d'),'.*','x','.^',int2str((i-1)));
	end
	[row,col]=size(str_f);
	str_f=strcat(str_f,')');
	f=inline(str_f,'x');
	cint=quadl(f,0,l_max,1e-10);
	%normalizing
	oa=oa/cint;
	str_f=strcat('x.*(x-',num2str(l_max),').*(');
	for i=1:(opt_n-1)
	str_f=strcat(str_f,num2str(oa(i),'%+d'),'.*','x','.^',int2str((i-1)));
	end
	[row,col]=size(str_f);
	str_f=strcat(str_f,')');
	f=inline(str_f,'x');
l=l_min:0.01:l_max;
poly_pdf_l=f(l);
figure
bar(C,N,'w')
hold on
plot(l,poly_pdf_l,'b--','LineWidth',2);
xlabel('Trace length','FontSize',15)
ylabel('迹长概率密度函数','FontSize',15)
legend('归一化统计直方图','最小二乘拟合密度函数',5)
legend boxoff
box off
hold off

[sme,poly_pdf_p]=Legendre_poly(opt_n,l_min,l_max,l,poly_pdf_l);
% %输出变量，每一个变量输出txt文件
poly_pdf_p%命名为“拟合多项式系数__从高阶到低阶”




%%%%%%%%%%%%%%%测窗迹线分维覆盖计算%%%%%%%%%%%%%%%%%%%%%
mwh=16;
mwv=8;
plane_p=load('测窗迹长端点坐标.txt'); %len_c_tra--6xn matrix, 迹长两端点坐标位于同一列数组中
%设定方格尺度向量
ms=[4 2 1 0.5 0.25];%方格尺度数组
figure
x=[0, 2*mwh, 2*mwh, 2*mwh, 2*mwh,     0,     0,  0];
y=[0,     0,     0, 2*mwv, 2*mwv, 2*mwv, 2*mwv,  0];
z=[0,     0,     0,    0,      0,     0,     0,  0];
colort=rand(1,3);
line(x,y,z,'Color',colort,'LineWidth',2.5)
[row col]=size(plane_p);
for i=1:col
    line(([plane_p(1,i),plane_p(4,i)]),([plane_p(2,i),plane_p(5,i)]),([plane_p(3,i),plane_p(6,i)]),'Color',colort,'LineWidth',2);
end
[rec_cv,rec_cc]=CalFractalDim2D(mwh,mwv,ms,plane_p);
axis equal
%输出txt文件
rec_cv %各个尺度下覆盖裂隙的方格数目，命名为“覆盖方格数”




% k_axis=2;%轴比
% beta_loc=45*pi/180;%裂隙中心到迹线垂线与长轴的夹角
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The following are the process to obtian PDF of the length of major axis %%%%%%%%%%%%%%%%%%%%%%%%%
% M=sqrt((cot(beta_loc))^2+1)/sqrt((k_axis*cot(beta_loc))^2+1);
% M=0.8771;
% %calculate the Abel integral equation
% a_max=l_max/M; %a_max--the maximum characteristic dimension
% a_min=5.5;%l_min/M;%a_min--the initial characteristic dimension
% syms a 
% miu_a_new=1;
% miu_a=0;
% % while ((abs(miu_a_new-miu_a)>0.5) && (a_min>0) )
% %     a_min=a_min-0.5
%     [miu_a,miu_a_new,g]=Joint_size_distri(opt_n,oa,a_min,l_max,M);
%     miu_a
%     miu_a_new
% % end
% g=vpa(g,4); 
% ai=a_min:0.05:a_max;
% [row,col]=size(ai);
% poly_pdf_ai=zeros(row,col);
% for i=1:col
% 	poly_pdf_ai(1,i)=subs(g,a,ai(i));
% end
% m=trapz(ai,poly_pdf_ai);
% g=g/m;
% % %输出变量为Matlab Data 文件
% save Size_PDF g a_max a_min
% 
% for i=1:col
% 	poly_pdf_ai(1,i)=subs(g,a,ai(i));
% end
% figure
% plot(ai,poly_pdf_ai,'k--','LineWidth',2,'MarkerSize',5)
% xlabel('椭圆裂隙长轴','FontSize',15)
% ylabel('概率密度','FontSize',15)
% legend('推断的尺寸分布函数',5)
% legend boxoff
% box off

