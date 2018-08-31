function [jg_dd_dip,es_k,flag]=EsFisherParaW(j_dd_dip,w)
%jg_dd_dip-2 x 1 matrix ，节理组的[倾向;倾角] 弧度
%j_dd_dip-2 x n matrix ，节理的[倾向;倾角] 弧度
%es_k- 估计的k值
%w-是权重
[row,col]=size(j_dd_dip);
l=0;m=0;n=0;
for j=1:col
	alpha=j_dd_dip(1,j); beta=j_dd_dip(2,j);
	%设置球极坐标，theta和phi
	theta=beta;
	if alpha<=0.5*pi
		phi=0.5*pi-alpha;
	else
		phi=2.5*pi-alpha;	
	end
	l=l+w(j)*sin(theta)*cos(phi);
	m=m+w(j)*sin(theta)*sin(phi);
	n=n+w(j)*cos(theta);
end
%合成距离
r=(l^2+m^2+n^2)^0.5;
%平均方向：mu是平均方向与z轴的夹角，nu是方位角
mu=acos(n/r); cos_phi=l/(r*sin(mu)); sin_phi=m/(r*sin(mu));
%cos_phi=-0.7071;sin_phi=-0.7071;
%根据象限判断，给出nu
if sin_phi>=0
	if cos_phi>=0
		nu=asin(sin_phi);	
	else
		nu=acos(cos_phi);
	end
else
	if cos_phi<=0
		nu=pi+asin(-sin_phi);	
	else
		nu=2*pi+asin(sin_phi);
	end
end
%mu=mu*180/pi; nu=nu*180/pi
%转换成地质上需要的倾向倾角,节理组joint group

if nu>0.5*pi
	m_dd=2.5*pi-nu;
else
	m_dd=0.5*pi-nu;
end
m_dip=mu;	
%给出节理组的产状
jg_dd_dip=[m_dd;m_dip];
%[row,col]=size(j_dd_dip);
%es_k=(col-1)/(col-r*col/sum(w));
%return

%%%%%以下报废

%转换矩阵nto
nto=[cos(mu)*cos(nu),cos(mu)*sin(nu),-sin(mu);-sin(nu),cos(nu),0;sin(mu)*cos(nu),sin(mu)*sin(nu),cos(mu);];
nto=PruneMartix(nto,1e-5);

c=zeros(1,col);
for j=1:col
	alpha=j_dd_dip(1,j); beta=j_dd_dip(2,j);
	%设置球极坐标，theta和phi，于青春教授的书
	theta=beta;
	if alpha<=0.5*pi
		phi=0.5*pi-alpha;
	else
		phi=2.5*pi-alpha;	
	end
	nv=nto*[(sin(theta)*cos(phi));(sin(theta)*sin(phi));cos(theta)];
	c(j)=1-nv(3,1);
end
%离散参数的估计
es_k0=col/sum(c);
es_k1=(col-1.0)*es_k0/col;
es_k2=(col-2.0)*es_k0/col;
es_k3=(1.0-1.0/col)^1.5*es_k0;
%统计检验是必要的卡方检验
%检验3个k的估计值，还是统计检验说的算~~~
%检验
%设定nk份，5份是较少的
nk=5; 
ndiv=(max(c)-min(c))/nk;
edge=zeros(1,(nk+2));
edge(1,1)=min(c);
for j=2:(nk+2)
	edge(1,j)=edge(1,1)+(j-1)*ndiv;
end
[nc,bin] = histc(c,edge); 

ts_k=[es_k1,es_k2,es_k3];
%记录chi值的矩阵
tp_k=zeros(1,3);
%求解理论计算概率 k*exp(-c*k);
p=zeros(1,(nk+1));
for i=1:3
	es_k=ts_k(1,i);
	for j=2:(nk+2)
		p(j-1)=-[exp(-es_k*edge(j))-exp(-es_k*edge(j-1))];
	end
	%计算卡方chi
	chi=0;
	for j=1:(nk+1)
		chi=chi+(nc(j)-col*p(j))^2/(col*p(j));
	end
	tp_k(1,i)=chi;
end
%找出最小的chi值
[chi,ind]=min(tp_k);
es_k=ts_k(1,ind);

%显著性水平alpha=0.05，chi2(1-alpha)=
if chi<chi2inv(0.95,(nk-2))
flag=1;
else
flag=0;
end