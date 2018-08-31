%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This function is used to calculate the polynomial form of trace %%%%%%
%%%% length using best square approchs , and the polynomial is enforced %%%
%%%% by Legendre orthogonal polynomial    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sme,poly_pdf_p]=Legendre_poly(n,l_min,l_max,l,pdf_l)

% n--the order of legendre polynomial n
% l_min,l_max---the miniomum and maximum trace length
% l--the discrete trace length
% pdf_l--the analytical PDF of trace length 

%%%%%%%%%%%%%%%%%%%%%%%% Calculating the polynomial form of trace length %%%%%%%%%%%%%%%%%%%%%%%%
a=zeros(1,(n+1));
dl=l(2)-l(1);
%linear transformate to interval [-1,1]
t=(l-0.5*(l_max+l_min))./(0.5*(l_max-l_min));
dt=t(2)-t(1);
[row,col]=size(t);
%P = legendre(n,X) computes the associated Legendre functions  of degree n
%and order m = 0,1,...,n, evaluated for each element of X. Argument n must be a scalar integer, and X must contain real values in the domain .
%If X is a vector, then P is an (n+1)-by-q matrix, where q = length(X). Each element P(m+1,i) corresponds to the associated Legendre function of degree n and order m evaluated at X(i).
for j=1:(n+1)
   p=legendre((j-1),t);
   f=p(1,:).*pdf_l;
   cint=(0.5*f(1)+sum(f(1,2:(col-1)))+0.5*f(col))*dt;
   a(j)=0.5*(2*j-1)*cint;  %it's the coefficient of the every Legendre orthogonal polynomial
end
%calculating the square root of error
f=pdf_l.*pdf_l;
%scalar product pdf with pdf
cint=(0.5*f(1)+sum(f(1,2:(col-1)))+0.5*f(col))*dt;
%root-mean-square error
sme=0;
for i=1:(n+1)
    sme=sme+2.0/(2*i-1)*a(i)^2;
end
sme=sqrt(cint-sme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficient of polynomial for the range of [-1 1]
cleg=zeros((n+1),(n+1));
for i=3:(n+1)
    cleg(1,(n+1))=1;
    cleg(2,n)=1;
    cleg(2,(n+1))=0;
    for j=1:(i-1)
        cleg(i,(n+1-i+j))=(2*i-3)/(i-1)*cleg((i-1),(n+1-(i-1)+j)); %Recursive relationship
    end
    for j=1:(i-2)
        cleg(i,(n+1-(i-2)+j))=cleg(i,(n+1-(i-2)+j))-(i-2)/(i-1)*cleg((i-2),(n+1-(i-2)+j));%Recursive relationship
    end
end
for i=1:(n+1)
    cleg(i,:)=a(i)*cleg(i,:);
end
clear p
p=sum(cleg); % coefficient of polynomial for the range of [-1 1]


%calculating coefficients through binomial theorem
clear a;  clear c;
cleg=zeros((n+1),(n+1));
a=2.0/(l_max-l_min);   b=(l_max+l_min)/(l_min-l_max);
for i=1:(n+1)
    k=n+1-i;
    for j=0:k
        cleg((n+2-i),(j+i))=p(i)*nchoosek(k,j)*(a^(k-j))*(b^j);
    end
end
%p are the coeffcients of the fitting ploynomial of the probability density function
%p(1) is the coeffcient of Highest order of the fitting ploynomial
p=sum(cleg);
f= polyval(p,l);
cint=(0.5*f(1)+sum(f(1,2:(col-1)))+0.5*f(col))*dl;
%poly_pdf_p--the fitted coefficient the polynomial, from high order to low order
poly_pdf_p=p/cint;%nomalization 

end