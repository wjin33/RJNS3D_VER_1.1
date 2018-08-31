function [insc_p,is_flag]=GetLSIP2D(line_p1,line_p2,poly_p1,poly_p2)
%line_p1,line_p2--3x1 matrix£¬the intersection line segment of the  circumscribed rectangular
%poly_p1,poly_p2--3x2 matrix£¬the line segment which consists one edge of the polygon
%insc_p--3x1 matrix£¬the coordinates of  intersection point
%is_flag--1:intersect£¬0:non-intersect

insc_p=[];
is_flag=lsLSIntersect2D(line_p1,line_p2,poly_p1,poly_p2);
if is_flag==1
u=((poly_p1(1,1)-line_p1(1,1))*(poly_p1(2,1)-poly_p2(2,1))-(poly_p1(2,1)-line_p1(2,1))*(poly_p1(1,1)-poly_p2(1,1)))/((line_p2(1,1)-line_p1(1,1))*(poly_p1(2,1)-poly_p2(2,1))-(line_p2(2,1)-line_p1(2,1))*(poly_p1(1,1)-poly_p2(1,1)));
insc_p(1,1)=line_p1(1,1)+u*(line_p2(1,1)-line_p1(1,1));
insc_p(2,1)=line_p1(2,1)+u*(line_p2(2,1)-line_p1(2,1));
insc_p(3,1)=0;
end
