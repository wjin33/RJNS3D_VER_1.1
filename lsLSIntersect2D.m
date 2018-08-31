function is_flag=lsLSIntersect2D(line_p1,line_p2,poly_p1,poly_p2)
%line_p1,line_p2--3x1 matrix£¬the intersection line segment of the  circumscribed rectangular
%poly_p1,poly_p2--3x2 matrix£¬the line segment which consists one edge of the polygon
%is_flag--1:intersection£¬0:non-intersect

t1=max(line_p1(1,1),line_p2(1,1));
t2=min(poly_p1(1,1),poly_p2(1,1));
t3=max(poly_p1(1,1),poly_p2(1,1));
t4=min(line_p1(1,1),line_p2(1,1));
t5=max(line_p1(2,1),line_p2(2,1));
t6=min(poly_p1(2,1),poly_p2(2,1));
t7=max(poly_p1(2,1),poly_p2(2,1));
t8=min(line_p1(2,1),line_p2(2,1));

mt1=((poly_p1(1,1)-line_p1(1,1))*(line_p2(2,1)-line_p1(2,1))-(line_p2(1,1)-line_p1(1,1))*(poly_p1(2,1)-line_p1(2,1)));
mt2=((line_p2(1,1)-line_p1(1,1))*(poly_p2(2,1)-line_p1(2,1))-(poly_p2(1,1)-line_p1(1,1))*(line_p2(2,1)-line_p1(2,1)));
mt3=((line_p1(1,1)-poly_p1(1,1))*(poly_p2(2,1)-poly_p1(2,1))-(poly_p2(1,1)-poly_p1(1,1))*(line_p1(2,1)-poly_p1(2,1)));
mt4=((poly_p2(1,1)-poly_p1(1,1))*(line_p2(2,1)-poly_p1(2,1))-(line_p2(1,1)-poly_p1(1,1))*(poly_p2(2,1)-poly_p1(2,1)));

%one line parallels with another, or superposition with each other's  part: non-intersect
if ((abs((t1-t4))<1e-5)&&(abs((t2-t3))<1e-5))|((abs((t6-t7))<1e-5)&&(abs((t5-t8))<1e-5))
	is_flag=0;
	return;
end

if (abs((t1-t2))<=1e-5)
t1=t2;
end
if (abs((t3-t4))<=1e-5)
t3=t4;
end
if (abs((t5-t6))<=1e-5)
t5=t6;
end
if (abs((t7-t8))<=1e-5)
t7=t8;
end

%decide the crossproduct is 0 or not,  tolerance is 1e-5
if (abs((mt1*mt2))<=1e-5)
mt1=0;
end
if (abs((mt3*mt4))<=1e-5)
mt3=0;
end

if (t1>=t2)&(t3>=t4)&(t5>=t6)&(t7>=t8)&((mt1*mt2)>=0)&((mt3*mt4)>=0)
is_flag=1;
else
is_flag=0;
end

end