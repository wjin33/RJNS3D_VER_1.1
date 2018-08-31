function j_dd_dip=GenFisherRand(jg_dd_dip,es_k,num_joint)
%jg_dd_dip--2x1, original trend and dip
%es_k--the estimated K value, pay attention to choose the appropriate k according to the number of sample
%num_joint--number the joints


%set the spherical coordinates parameter, theta and phi
mu=jg_dd_dip(2,1);
if jg_dd_dip(1,1)<=0.5*pi
	nu=0.5*pi-jg_dd_dip(1,1);
else
	nu=2.5*pi-jg_dd_dip(1,1);
end

nto=[cos(mu)*cos(nu),cos(mu)*sin(nu),-sin(mu);
     -sin(nu),cos(nu),0;
     sin(mu)*cos(nu),sin(mu)*sin(nu),cos(mu);];
otn=inv(nto); otn=PruneMartix(otn,1e-5);
rand('state',sum(100*clock));
j_dd_dip=zeros(2,num_joint);
c=zeros(1,num_joint);
for j=1:num_joint
	%generate rand fisher parameter in spherical coordinates
	theta_t=acos((log((exp(es_k)-(exp(es_k)-exp(-es_k))*rand))/es_k));
	phi_t=2*pi*rand; 	c(j)=1-cos(theta_t);
	%transform to N-Y,E-X globle coordinates
	l=otn(1,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
	m=otn(2,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
	n=otn(3,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
	theta=acos(n);
	while theta>0.5*pi
		theta_t=acos((log((exp(es_k)-(exp(es_k)-exp(-es_k))*rand))/es_k));
		phi_t=2*pi*rand; 	c(j)=1-cos(theta_t);
		%transform to N-Y,E-X globle coordinates
		l=otn(1,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
		m=otn(2,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
		n=otn(3,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
		l=PruneMartix(l,1e-5); m=PruneMartix(m,1e-5); n=PruneMartix(n,1e-5);
		theta=acos(n);
    end
    dip=theta;
    
    
	cos_phi=l/sin(theta); sin_phi=m/sin(theta); 
	%according to quandrant, choose phi
	if sin_phi>=0
		if cos_phi>=0
			phi=asin(sin_phi);
		else
			phi=acos(cos_phi);
		end
	else
		if cos_phi<=0
			phi=pi+asin(-sin_phi);	
		else
			phi=2*pi+asin(sin_phi);
		end
	end
	%transform to trend and dip
	if phi>0.5*pi
		dd=2.5*pi-phi;
	else
		dd=0.5*pi-phi;
	end
		
	j_dd_dip(1,j)=dd;
	j_dd_dip(2,j)=dip;
end

%the chi-square statistical test
%at least 5 intervals
nk=5; 
ndiv=(max(c)-min(c))/nk;
edge=zeros(1,(nk+1));
edge(1,1)=min(c);
for j=2:(nk+1)
	edge(1,j)=edge(1,1)+(j-1)*ndiv;
end
[nc, bin] = histc(c,edge); 
%calculater the theoretical probability:      k*exp(-c*k);
p=zeros(1,(nk+1));
for j=2:(nk+1)
	p(j-1)=-[exp(-es_k*edge(j))-exp(-es_k*edge(j-1))];
end
%calculate chi
chi=0;
for j=1:nk
	chi=chi+(nc(j)-num_joint*p(j))^2/(num_joint*p(j));
end
%significance level alpha=0.05
if chi<chi2inv(0.95,(nk-1))
    flag=1;
else
    flag=0;
end


%if the rand number cann't satisfy chi-square statistical test, then do it again
while flag==0
	rand('state',sum(100*clock));
	j_dd_dip=zeros(2,num_joint);
	c=zeros(1,num_joint);
	for j=1:num_joint
		theta_t=acos((log((exp(es_k)-(exp(es_k)-exp(-es_k))*rand))/es_k));
		phi_t=2*pi*rand; 	c(j)=1-cos(theta_t);
		l=otn(1,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
		m=otn(2,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
		n=otn(3,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
		theta=acos(n);
		while theta>0.5*pi
			theta_t=acos((log((exp(es_k)-(exp(es_k)-exp(-es_k))*rand))/es_k));
			phi_t=2*pi*rand; 	c(j)=1-cos(theta_t);
			l=otn(1,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
			m=otn(2,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
			n=otn(3,:)*[sin(theta_t)*cos(phi_t);sin(theta_t)*sin(phi_t);cos(theta_t)];
			l=PruneMartix(l,1e-5); m=PruneMartix(m,1e-5); n=PruneMartix(n,1e-5);
			theta=acos(n);
		end
		cos_phi=l/sin(theta); sin_phi=m/sin(theta); 
		if sin_phi>=0
			if cos_phi>=0
				phi=asin(sin_phi);	
			else
				phi=acos(cos_phi);
			end
		else
			if cos_phi<=0
				phi=pi+asin(-sin_phi);	
			else
				phi=2*pi+asin(sin_phi);
			end
        end
		if phi>0.5*pi
			dd=2.5*pi-phi;
		else
			dd=0.5*pi-phi;
		end
		dip=theta;	
		j_dd_dip(1,j)=dd;
		j_dd_dip(2,j)=dip;
    end
    
    %the chi-square statistical test
    %at least 5 intervals
    nk=5; 
    ndiv=(max(c)-min(c))/nk;
    edge=zeros(1,(nk+1));
    edge(1,1)=min(c);
    for j=2:(nk+1)
        edge(1,j)=edge(1,1)+(j-1)*ndiv;
    end
    [nc, bin] = histc(c,edge); 
    %calculater the theoretical probability:      k*exp(-c*k);
    p=zeros(1,(nk+1));
    for j=2:(nk+1)
    	p(j-1)=-[exp(-es_k*edge(j))-exp(-es_k*edge(j-1))];
    end
    %calculate chi
    chi=0;
    for j=1:nk
    	chi=chi+(nc(j)-num_joint*p(j))^2/(num_joint*p(j));
    end
    %significance level alpha=0.05
    if chi<chi2inv(0.95,(nk-1))
      flag=1;
    else
       flag=0;
    end
end
end