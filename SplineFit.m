clear all
clc
%theta--the angle bewteen trace and horizontal edge of the measure window
%set len_c_tra: length of trace inside the MW;
len_c_tra=[0.57428	0.23345	0.7912	0.24166	0.99247	0.59506	0.68264	0.95525	0.32311	0.49244	0.31385	0.48052	0.49244	0.70235	0.6261	1.9721 0.30265	0.52	0.35847	0.08544	0.99045	0.38013	0.45255	0.4993	0.60033	0.536	0.20616	0.20616	0.59363	0.55946	0.46098	1.321	0.98615	0.74027	1.1314	0.77026	1.1536;];
theta=1.1699;
mwh=5;
mwv=1;
es_dens=2.7418;% density of trace in MW
len_c_tra = len_c_tra(:);
[F_,X_] = ecdf(len_c_tra,'Function','cdf');  % compute empirical cumulative distribution function cdf
%[f,x] = ecdf(y) calculates the Kaplan-Meier estimate of the cumulative distribution function (cdf),
%also known as the empirical cdf. y is a vector of data values. f is a vector of values of the empirical cdf evaluated at x.
Bin_.rule = 3;
Bin_.nbins = 3;
[C_,E_] = dfswitchyard('dfhistbins',len_c_tra,[],[],Bin_,F_,X_);%switchyard for Distribution Fitting.
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
%n = ecdfhist(f, x) takes a vector f of empirical cumulative distribution function (cdf) values and a vector x of evaluation points,
%and returns a vector n containing the heights of histogram bars for 10 equally spaced bins. The function computes the bar heights 
%from the increases in the empirical cdf, and normalizes them so that the area of the histogram is equal to 1. 
tyc=[]; 
tdh=[];
yc=C_;
nc=N_*length(len_c_tra)/sum(N_);
dh=((nc./es_dens)./(cos(theta)*sin(theta).*(yc.^2)-(2*mwh*sin(theta)+2*mwv*cos(theta)).*yc+4*mwh*mwv))./(yc(2)-yc(1));
tyc=[tyc,yc]; tdh=[tdh,dh];
dmax=max(len_c_tra)+0.01;
% dmax=2*mwv/sin(theta)
% dmax=2.06
%try to use spline to fit
tyc=[tyc,yc]; tdh=[tdh,dh];
[bm1,ix] = sort(tyc(1,:),'ascend');
tyc=tyc(1,ix);  tdh=tdh(1,ix);
tyc=[0,tyc,dmax]; tdh=[0,tdh,0]
%Spline of PDF
s=csapi(tyc,tdh);

x=0:0.001:dmax;
y=fnval(s,x);
cint=trapz(x,y)
s.coefs=s.coefs/cint
%verifying by caiculus
y=fnval(s,x);
cint=trapz(x,y)
x=0:0.002:dmax;
y=fnval(s,x);
%probability distribution function of trace
trmax=dmax; trmin=0;
l=trmin:0.002:trmax;
[row,col]=size(l);
pdf=zeros(row,col);
for j=1:col
    pdf(j)=fnval(s,l(j));
end
hold on
plot(l,pdf,'m-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the order of legendre polynomial n
n=16;
a=zeros(1,(n+1)); 
x=0:0.002:trmax;
%linear transformate to interval [-1,1]
t=(x-0.5*(trmax+trmin))./(0.5*(trmax-trmin));
dt=t(2)-t(1);
%P = legendre(n,X) computes the associated Legendre functions  of degree n
%and order m = 0,1,...,n, evaluated for each element of X. Argument n must be a scalar integer, and X must contain real values in the domain . If X is a vector, then P is an (n+1)-by-q matrix, where q = length(X). Each element P(m+1,i) corresponds to the associated Legendre function of degree n and order m evaluated at X(i).
for j=1:(n+1)
   p=legendre((j-1),t);
   f=p(1,:).*pdf;    %Best Square Approximation
   inttz=f(1)+f(col);
   for i=2:(col-1)
      inttz=inttz+2*f(i);
   end
   inttz=0.5*inttz*dt;
   a(j)=0.5*(2*j-1)*inttz;
end


[row,col]=size(t);
f=pdf.*pdf;
%scalar product pdf with pdf
inttz=f(1)+f(col);
for i=2:(col-1)
    inttz=inttz+2*f(i);
end
inttz=0.5*inttz*dt;
%root-mean-square error
sme=0;
for i=1:(n+1)
    sme=sme+2.0/(2*i-1)*a(i)^2;
end
sme=sqrt(inttz-sme)
f=0;
for j=1:(n+1)
    p=legendre((j-1),t);
    f=f+a(j)*p(1,:);
end



%% coefficient of Legendre orthogonal polynomial
cleg=zeros((n+1),(n+1));
for i=3:(n+1)
    cleg(1,(n+1))=1;
    cleg(2,n)=1;
    cleg(2,(n+1))=0;
    for j=1:(i-1)
        cleg(i,(n+1-i+j))=(2*i-3)/(i-1)*cleg((i-1),(n+1-(i-1)+j));
    end
    for j=1:(i-2)
        cleg(i,(n+1-(i-2)+j))=cleg(i,(n+1-(i-2)+j))-(i-2)/(i-1)*cleg((i-2),(n+1-(i-2)+j));
    end
end

for i=1:(n+1)
    cleg(i,:)=a(i)*cleg(i,:);
end



clear p
p=sum(cleg);
%verify the polynomial 
f=0;
for j=1:(n+1)
    f=f+p(j)*t.^(n+1-j);
end

[row,col]=size(t);
inttz=f(1)+f(col);
for i=2:(col-1)
    inttz=inttz+2*f(i);
end
inttz=0.5*(trmax-trmin)*0.5*inttz*dt
p=p/inttz;%being  normalized

%acquire coefficient through binomial theorem
clear a;  clear c;
cleg=zeros((n+1),(n+1));
a=2.0/(trmax-trmin);   b=(trmax+trmin)/(trmin-trmax);
for i=1:(n+1)
    k=n+1-i;
    for j=0:k
        cleg((n+2-i),(j+i))=p(i)*nchoosek(k,j)*(a^(k-j))*(b^j);
    end
end
p=sum(cleg);
%p are the coeffcients of the fitting ploynomial of the probability density function
%p(1) is the coeffcient of Highest order of the fitting ploynomial
[row,col]=size(p);
p(row,col)=0;
f=polyval(p,x);
[row,col]=size(t);
inttz=f(1)+f(col);
for i=2:(col-1)
    inttz=inttz+2*f(i);
end
inttz=0.5*(trmax-trmin)*0.5*inttz*dt
p=p/inttz;
x=0:0.002:trmax;
f=polyval(p,x);


plot(x,f,'r--','LineWidth',3)
hold on
plot(x,pdf,'b-','LineWidth',2)
xlabel('Trace Length','FontSize',15)
ylabel('Probability Density Function','FontSize',15)
set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9],'XTickLabel',{'0','1','2','3','4','5','6','7','8','9'},'FontSize',10)
legend('Best polynomial approximation of PDF','Spline of distribution of trace',2)
legend boxoff
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The following are the process to obtian PDF of radius%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the Abel integral equation
a=p;
[row,col]=size(a);
ind=(col-1):-1:1;
a=a(1,ind);

syms x y F
F=sym('0');
for k=1:(col-1)
	f=a(k)*(k-1)*(y^(k-2))*(y^2-x^2)^(-1/2);
	F=F+int(f,y,x,trmax);
end
F=-2*F*x/pi;


ix=[eps:0.002:trmax];
[row,col]=size(ix);
iy=zeros(row,col);
for i=1:col
	iy(1,i)=subs(F,x,ix(i));
end
%Z = trapz(X,Y) computes the integral of Y with respect to X using trapezoidal integration. 
m=1.0/trapz(ix,iy)
hold on
%distribution of radius; <<<<refer to P72 _doctoral dissertation>>>>>>
F=m*F;
for i=1:col
	iy(1,i)=ix(i)*subs(F,x,ix(i));
end
mc=trapz(ix,iy)

%calculateing PDF of radius
for i=1:col
	iy(1,i)=subs(F,x,ix(i));
end
figure(2)
plot(ix,iy,'ro--','LineWidth',2,'MarkerSize',5)
hold on
plot(ix,0.125*ones(row,col),'b-','LineWidth',2,'MarkerSize',5)
xlabel('Radius','FontSize',15)
ylabel('PDF of radius','FontSize',15)
axis([-1,9,0,0.2])
set(gca,'XTick',[-1,0,1,2,3,4,5,6,7,8,9],'XTickLabel',{'-1','0','1','2','3','4','5','6','7','8','9'},'FontSize',20)
legend('Inferential value','Theoretical value',5)
legend boxoff
box off