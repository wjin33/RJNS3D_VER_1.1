function on_flag=IsPointOnLS2D(p,poly_p1,poly_p2) 
%poly_p1,poly_p2--3x2 matrix£¬the line segment which is one edge of the polygon
%p--3x1 matrix, coordinate of the point
%on_flag--if it equals 1:within the line;  0: outside the line
on_flag=0;
%controling the coordinate of P to be in the rectangular that Diagonal vertex is P1, P2; Including the speacial cases of horizontal and vertical
if ((p(1,1)>min(poly_p2(1,1),poly_p1(1,1)))|(abs((p(1,1)-min(poly_p2(1,1),poly_p1(1,1))))<=1e-5))&...
((p(1,1)<max(poly_p2(1,1),poly_p1(1,1)))|(abs((p(1,1)-max(poly_p2(1,1),poly_p1(1,1))))<=1e-5))&...
((p(2,1)>min(poly_p2(2,1),poly_p1(2,1)))|(abs((p(2,1)-min(poly_p2(2,1),poly_p1(2,1))))<=1e-5))&...
((p(2,1)<max(poly_p2(2,1),poly_p1(2,1)))|(abs((p(2,1)-max(poly_p2(2,1),poly_p1(2,1))))<=1e-5))
	cross_flag=(poly_p2(1,1)-poly_p1(1,1))*(p(2,1)-poly_p1(2,1))-(poly_p2(2,1)-poly_p1(2,1))*(p(1,1)-poly_p1(1,1));  %cross product	`
	if abs(cross_flag)<=1e-3      %the tolerance 
		on_flag=1;
	else
		on_flag=0;
	end
end