function is_loc_flag=LSTraLoc2D(nvm,line_p1,line_p2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%precondition;the trace must intersect with the ploygon before use this function%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line_p1,line_p2--3x1 matrix£¬the intersection line segment, 2-dimension
%nvm--3xn matri, the vertexes of polygon, stored in columns;and the polygon will be divided
%is_loc_flag--0:inside the polygon, 1:line segment is cut, 2:line segment cross the ploygon
%t_m--3xn matrix,store the coordinates temporarily 

nvm=[nvm,nvm(:,1)];
[row,col]=size(nvm);
t_m=[];
for i=1:(col-1)
   poly_p1=nvm(:,i);   poly_p2=nvm(:,(i+1));
   [insc_p,is_flag]=GetLSIP2D(line_p1,line_p2,poly_p1,poly_p2);
	if is_flag==1
		t_m=[t_m,insc_p];
	end
end
[row,col]=size(t_m);
if col<=1
	is_loc_flag=col;
%find if the intersected point is located on both endpoint of two adjacent lines of polygon
elseif ((t_m(1,1)-t_m(1,2))^2+( t_m(2,1)-t_m(2,2))^2)^0.5<1e-5  
	is_loc_flag=1;
else
	is_loc_flag=2;
end