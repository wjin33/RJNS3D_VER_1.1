function point_c=PoissonProcess3D2(lumtav,dx,dy,dz)
%lumtav--the intensity of Poisson points 3 dimension
%lumtal--the intensity of Poisson points 1 dimension
%dx--length  dy--width(north direction)  dz--high
%pay attention to the unit
%point_c--3xn matrix, the coordinate of center of the joint
%div_len--1x1 matrix, the number the simulated joint

rand('state',sum(100*clock))%generate the exponent random number
lumtal=lumtav*dx*dz;%the line density of Poisson points alone y-axis,lumtal=1.0/traex; Maximum Likelihood
traex=1.0/lumtal;%the number of Poisson points in unit
%generate the showing section length
div_len=round((lumtal*dy));
flag=1;
while flag
    %exprnd(MU,m,n) generates exponential random numbers with mean MU, where scalars m and n are the row and column dimensions of R.
    trnd=exprnd(1.0/lumtal,1,div_len);  
    slen=sum(trnd);
    if (slen<dy)
        flag=0;
    else
        flag=1;
    end
end
point_c=zeros(3,div_len);
%pay attention to the axis cut by  measure window, and main showing direction
for i=1:div_len
    point_c(1,i)=dx*rand;
    point_c(2,i)=sum(trnd(1,1:i));
    point_c(3,i)=dz*rand;
end
