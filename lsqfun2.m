function [oa,res,aic,f]=lsqfun2(n,ix,iy,d)
%n--the oder the function
%ix--1xn;  iy--1xn
%oa--the coefficient the polynomial, from low to high
%str_f--the string of the polynomial
%res--residuals
%aic--smaller, better
%d--diameter
str_f=strcat('x.*(x-',num2str(d),').*(');
for i=1:n-1
str_f=strcat(str_f,'a(',int2str(i),')','.*','x','.^',int2str(i-1),'+');
end
[row,col]=size(str_f);
str_f(1,col)=')'; %replace the last plus sign(+) with ')'
f=inline(str_f,'a','x');
ini_a=1.0*ones(1,(n-1)); %initial fitting parameter
%beta = nlinfit(X,y,fun,beta0) estimates the coefficients of a nonlinear function using least squares. 
%y is a vector of response (dependent variable) values. Typically, X is a design matrix of predictor (independent variable) values, 
%with one row for each value in y. However, X can be any array that fun can accept. fun is a function of the form yhat = myfun(beta,X)
%where beta is a coefficient vector, and X is the design matrix. fun returns a vector yhat of fitted y values. 
%beta0 is a vector containing initial values for the coefficients. [beta,r,J] = nlinfit(X,y,fun,beta0) returns the fitted coefficients, beta, 
%the residuals, r, and the Jacobian, J.
[oa,res]=nlinfit(ix,iy,f,ini_a);
res=sum((res.*res));
[row,col]=size(ix);
sigma2=res/col;
%%what's this??????
aic=col*(log10((2*pi*sigma2))+1)+2*n;
end