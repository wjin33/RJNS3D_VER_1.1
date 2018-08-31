function grand=JointSizeRand_TraFractal(D,a_min,a_max,numr)
%a_max--the maximum major axis 
%a_min--the minimum major axis 
%grand--
%numr--the number of major axis need be generated

syms x;
g=(D+1)*a_min^(D+1)*x^(-D-2);

%calculate the maximum value
ix=[a_min:0.001:a_max];
[row,col]=size(ix);
iy=zeros(row,col);
for i=1:col
	iy(1,i)=subs(g,x,ix(i));
end
c=max(iy);

%generate random number between a_min and a_max
rand('state',sum(100*clock))
y=a_min+(a_max-a_min)*rand;
u=rand;
gy=subs(g,x,y);
grand=zeros(1,numr);
flag=1;
while flag<=numr
	if u<=(gy/c)
		grand(1,flag)=y;
		flag=flag+1;
	end
	y=a_min+(a_max-a_min)*rand;
	gy=subs(g,x,y);
	u=rand;
end

end