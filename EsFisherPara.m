function [jg_dd_dip,es_k,flag]=EsFisherPara(j_dd_dip)

%j_dd_dip--2xn matrix, trend and dip of the investigated joint
%jg_dd_dip--2x1, the mean trend and dip
%es_k--the estimated K value,



[row,col]=size(j_dd_dip);
l=0;m=0;n=0;

%%

for j=1:col
	alpha=j_dd_dip(1,j); beta=j_dd_dip(2,j);
	%set the spherical coordinates parameter, theta and phi
	theta=beta;
	if alpha<=0.5*pi
		phi=0.5*pi-alpha;
	else
		phi=2.5*pi-alpha;	
	end
	l=l+sin(theta)*cos(phi);
	m=m+sin(theta)*sin(phi);
	n=n+cos(theta);
end

r=(l^2+m^2+n^2)^0.5; %length of the mean direction
%mean direction: mu is the angle between z-axis and r, phi is the  the angle between x-axis and the projection of r in xy plane
mu=acos(n/r); cos_phi=l/(r*sin(mu)); sin_phi=m/(r*sin(mu));

%according to quandrant, choose phi
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

%transform to trend and dip for joint group
if nu>0.5*pi
	m_dd=2.5*pi-nu;
else
	m_dd=0.5*pi-nu;
end
m_dip=mu;	
jg_dd_dip=[m_dd;m_dip];
%%


nto=[cos(mu)*cos(nu),cos(mu)*sin(nu),-sin(mu);
     -sin(nu),cos(nu),0;
     sin(mu)*cos(nu),sin(mu)*sin(nu),cos(mu);];
nto=PruneMartix(nto,1e-5);

c=zeros(1,col);
for j=1:col
	alpha=j_dd_dip(1,j); beta=j_dd_dip(2,j);
	%%set the spherical coordinates parameter, theta and phi, refer to the
	%%book by Prof.Yu Qingchun
	theta=beta;
	if alpha<=0.5*pi
		phi=0.5*pi-alpha;
	else
		phi=2.5*pi-alpha;	
	end
	nv=nto*[(sin(theta)*cos(phi));(sin(theta)*sin(phi));cos(theta)];
	c(j)=1-nv(3,1);
end
%estimate the discrete parameter
es_k0=col/sum(c);
es_k1=(col-1.0)*es_k0/col;
es_k2=(col-2.0)*es_k0/col;
es_k3=(1.0-1.0/col)^1.5*es_k0;
%the chi-square statistical test 
%at least 5 intervals
nk=5; 
ndiv=(max(c)-min(c))/nk;
edge=zeros(1,(nk+1));
edge(1,1)=min(c);
for j=2:(nk+1)
	edge(1,j)=edge(1,1)+(j-1)*ndiv;
end
[nc,bin] = histc(c,edge); 

ts_k=[es_k1,es_k2,es_k3];
tp_k=zeros(1,3); %store chi for each es_k
p=zeros(1,(nk+1));
for i=1:3       
	es_k=ts_k(1,i);
	for j=2:(nk+1)  %calculater the theoretical probability:      k*exp(-c*k);
		p(j-1)=-[exp(-es_k*edge(j))-exp(-es_k*edge(j-1))];
    end
	chi=0;
	for j=1:nk %calculate chi for each es_k
		chi=chi+(nc(j)-col*p(j))^2/(col*p(j));
	end
	tp_k(1,i)=chi;
end
%choose the minimal es_k
[chi,ind]=min(tp_k);
es_k=ts_k(1,ind);

%significance level alpha=0.05
if chi<chi2inv(0.95,(nk-1))
    flag=1;
else
    flag=0;
end

end