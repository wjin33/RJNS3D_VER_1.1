function [rec_p,lp_flag]=GetPolyLSIP2D(line_p1,line_p2,nvm)
%line_p1,line_p2--3x1 matrix£¬the intersection line segment
%nvm--3xn matri, the vertexes of polygon, stored in columns
%insc_p--3x1 matrix£¬the coordinates of  intersection point
%lp_flag--0:inside the polygon, 1:line segment is cut, 2:line segment cross the ploygon, -1:outside the polygon
%rec_p--6x1 matrix -coordinates of the the coordinates of  intersection point and the pints inside the polygon

rec_p=[];
nvm=[nvm, nvm(:,1)];
[row,col]=size(nvm);
lp_flag=0; 
%calculate the position of endpoints of the line
in_flag1=IsPointInPoly2D(nvm,line_p1); in_flag2=IsPointInPoly2D(nvm,line_p2); 

%inside the polygon
if ((in_flag1==1)||(in_flag1==0))&&((in_flag2==1)||(in_flag2==0))
	lp_flag=0;
	return
end

%have one intersection point between the polygon and the line
if ((in_flag1==1)&&(in_flag2==-1))||((in_flag1==-1)&&(in_flag2==1))
	rec_p=[];
	for i=1:(col-1)
		poly_p1=nvm(:,i);
		poly_p2 =nvm(:,(i+1));
		[insc_p,is_flag]=GetLSIP2D(line_p1,line_p2,poly_p1,poly_p2);
		if (is_flag==1)  
			lp_flag=1;
			if in_flag1==1	
				rec_p=[insc_p;line_p1];
			else
				rec_p=[insc_p;line_p2];
			end
			return;
		end
	end
end

%one of the endpoints is on the polygon
if ((in_flag1==0)&&(in_flag2==-1))||((in_flag1==-1)&&(in_flag2==0))
	count=0;	rec_p=[];
	for i=1:(col-1)
		poly_p1=nvm(:,i);
		poly_p2 =nvm(:,(i+1));
		[insc_p,is_flag]=GetLSIP2D(line_p1,line_p2,poly_p1,poly_p2);
		if (is_flag==1)  
			count=count+1;
			rec_p=[rec_p;insc_p];
		end
	end

	if count==1
		lp_flag=-1;
		return;
	else
		lp_flag=2;
		if ((rec_p(1,1)-rec_p(4,1))^2+(rec_p(2,1)-rec_p(5,1))^2)^0.5<1e-5
			lp_flag=1;
			return;
		end 
		return;
	end	
end

%both endpoints are outside the polygon
if (in_flag1==-1)&&(in_flag2==-1)
	count=0;  rec_p=[];
	for i=1:(col-1)
		poly_p1=nvm(:,i);
		poly_p2 =nvm(:,(i+1));
		[insc_p,is_flag]=GetLSIP2D(line_p1,line_p2,poly_p1,poly_p2);
		if (is_flag==1)  
			count=count+1;
			rec_p=[rec_p;insc_p];
		end
	end
if count==0
	lp_flag=-1;  %outside the polygon
	return;
else
	lp_flag=2;   %cross the polygon, but one of the endpoints is on the polygon
	return;
end

end