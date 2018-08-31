%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This function is used to calculate the polynomial form of joint size(major axis) %%%%%%
%%%% using the proposed approchs in the artical             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [miu_a,miu_a_new,g]=Joint_size_distri(n,poly_pdf_lc,a_min,l_max,M)

% n--the order of legendre polynomial n
% a_min--the miniomum  length of major axis
% l_max---the miniomum and maximum trace length
% l--the discrete trace length
% pdf_l--the analytical PDF of trace length 

%calculating the analytical formula of characteristic dimension for elliptical joints
%string of the  polynomial function
syms L ;
cs=strcat('L*(L-',num2str(l_max),')*(');
for i=1:(n-1)
	cs=strcat(cs,num2str(poly_pdf_lc(i),'+%d'),'*','L','^',int2str((i-1)));
end
cs=strcat(cs,')');
cs=strcat('pf=',cs);
eval(cs);  %EVAL(s), where s is a string, causes MATLAB to execute the string as an expression or statement.
coef=sym2poly(pf);   % SYM2POLY(P) returns a row vector containing the coefficients of the symbolic polynomial.
%obtain the coefficients in up power arrangement
[row,col]=size(coef);
ind=(col-1):-1:1;
coef=coef(1,ind);
%g--elliptical joint size distribution function, it's a symbolic function
a_max=l_max/M;
syms a G
G=sym('0');  %x = sym('x') creates the symbolic variable with name 'x' and stores the result in x.
for k=1:(col-1)
	f=coef(k)*(k-1)*(L^(k-2))*(L^2-(M*a)^2)^(-1/2);
	G=G+int(f,L,M*a,l_max);  
end
g0=-2*G*M*a/pi;%expression lack the first moment
g0=vpa(g0,4);
ia=a_min:0.05:a_max;
[row,col]=size(ia);
iy=zeros(row,col);
for i=1:col
	iy(1,i)=subs(g0,a,ia(i));
end
%using Trapezoidal method to calculate the integral
dh=ia(2)-ia(1);
cint=(0.5*iy(1)+sum(iy(1,2:(col-1)))+0.5*iy(col))*dh;
%First order moment of size of elliptical joints
miu_a=1/cint;
g=miu_a*g0;
g1=miu_a*g0*a;
g1=vpa(g1,4);
for i=1:col
	iy(1,i)=subs(g1,a,ia(i));
end
%using Trapezoidal method to calculate the integral
miu_a_new=(0.5*iy(1)+sum(iy(1,2:(col-1)))+0.5*iy(col))*dh;

end