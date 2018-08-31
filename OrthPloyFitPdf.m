%according to the theoretical distribution by  S. D. Priest(2004)
trmax=8; trmin=0;% maximum and minimum trace
l=0:0.002:trmax;
l(1,1)=0.000001;
[row,col]=size(l);
pdf=zeros(row,col);
for j=1:col
    pdf(j)=2*l(j)*(log(trmax+sqrt(trmax^2-l(j)^2))-log(l(j)))/(trmax^2-trmin^2);
end
hold on
plot(l,pdf,'m-')

%    dx=0.002;
%    inttz=pdf(1)+pdf(col);
%    for i=2:(col-1)
%       inttz=inttz+2*pdf(i);
%    end
%    inttz=0.5*inttz*dx

   
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
   f=p(1,:).*pdf;   
   inttz=f(1)+f(col);
   for i=2:(col-1)
      inttz=inttz+2*f(i);
   end
   inttz=0.5*inttz*dt;
   a(j)=0.5*(2*j-1)*inttz;  %it's the coefficient of the every Legendre  polynomial
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating the square root of error
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
% f=0;
% for j=1:(n+1)
%     p=legendre((j-1),t);
%     f=f+a(j)*p(1,:);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
% f=0;
% for j=1:(n+1)
%     f=f+p(j)*t.^(n+1-j);
% end
% [row,col]=size(t);
% inttz=f(1)+f(col);
% for i=2:(col-1)
%     inttz=inttz+2*f(i);
% end
% inttz=0.5*(trmax-trmin)*0.5*inttz*dt
% p=p/inttz;%being  normalized

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
% p(row,col)=0;
% f=polyval(p,x);
% [row,col]=size(x);
% dx=x(2)-x(1);
% inttz=f(1)+f(col);
% for i=2:(col-1)
%     inttz=inttz+2*f(i);
% end
% inttz=0.5*inttz*dx
% p=p/inttz;
% x=0:0.002:trmax;
% f=polyval(p,x);


plot(x,f,'r--','LineWidth',3)
hold on
plot(x,pdf,'b-','LineWidth',2)
xlabel('Trace Length','FontSize',15)
ylabel('Probability Density Function','FontSize',15)
set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9],'XTickLabel',{'0','1','2','3','4','5','6','7','8','9'},'FontSize',10)
legend('Best polynomial approximation of PDF','Theoretical value of distribution of trace',2)
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
