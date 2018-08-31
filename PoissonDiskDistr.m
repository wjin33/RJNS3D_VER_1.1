function [m,mc,g]=PoissonDiskDistr(oa,pfn,dmax,dstart)
%(pfn+1)--the oder of the function which need to be fitted
%g--diameter distribution function, it's a symbolic function
%dsmax--the maximum diameter 
%dstart--the initial diameter
%oa--the fitted coefficient the polynomial, from low to high
%m--A order moments
%mc--A order moments

%string of the  polynomial function
syms y ;
cs=strcat('y*(y-',num2str(dmax),')*(');
for i=1:(pfn-1)
	cs=strcat(cs,num2str(oa(i),'+%d'),'*','y','^',int2str((i-1)));
end
[row,col]=size(cs);
cs=strcat(cs,')');
cs=strcat('pf=',cs);
eval(cs);  %EVAL(s), where s is a string, causes MATLAB to execute the string as an expression or statement.
a=sym2poly(pf);    % SYM2POLY(P) returns a row vector containing the coefficients of the symbolic polynomial P.
%botain the coefficients in up power arrangement
[row,col]=size(a);
ind=(col-1):-1:1;
a=a(1,ind);
%calculate the Warburton formula
syms x g0
g0=sym('0');  %x = sym('x') creates the symbolic variable with name 'x' and stores the result in x.
for k=1:(col-1)
	f=a(k)*(k-1)*(y^(k-2))*(y^2-x^2)^(-1/2);
	g0=g0+int(f,y,x,dmax);
end
%expression lack the first moment
g0=-2*g0*x/pi;


ix=[dstart:0.001:dmax];
[row,col]=size(ix);
iy=zeros(row,col);
for i=1:col
	iy(1,i)=subs(g0,x,ix(i));
end
%using Trapezoidal method to calculate the integral
dh=ix(2)-ix(1);
cint=(0.5*iy(1)+sum(iy(1,2:(col-1)))+0.5*iy(col))*dh;
m=1/cint;%A order moments

%refer to the essay
g=m*g0;
g1=m*g0*x;
for i=1:col
	iy(1,i)=subs(g1,x,ix(i));
end
%using Trapezoidal method to calculate the integral
mc=(0.5*iy(1)+sum(iy(1,2:(col-1)))+0.5*iy(col))*dh;
end