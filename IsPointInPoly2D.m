function in_flag=IsPointInPoly2D(nvm,p) 
%nvm--the vertexes of polygon, stored in columns
%p--3x1 matrix£¬the coordinate of P
%in_flag--equal 1:inside the polygon;  equal 0:on the polygon;  equal -1:outside the polygon
%poly_p1,poly_p2--3x2 matrix£¬the line segment which consists one edge of the polygon
[row,col]=size(nvm);
nvm=nvm(:,col:-1:1);
nvm=[nvm, nvm(:,1)];
[row,col]=size(nvm);
count = 0; 
inf_p=zeros(3,1);
inf_p(2,1)=p(2,1);
inf_p(1,1)=-1e10;
for i=1:(col-1)
	poly_p1=nvm(:,i);
	poly_p2 =nvm(:,(i+1));
	if (IsPointOnLS2D(p, poly_p1,poly_p2)==1)  
    	in_flag=0;
	return;
	end
	if (abs((poly_p1(2,1)-poly_p2(2,1))) < 1e-3 ) %horizontal line
	continue; 
	end
	if (IsPointOnLS2D(poly_p1,p,inf_p)==1) 
		if( poly_p1(2,1) > poly_p2(2,1) ) 
		count=count+1; 
		end
	elseif (IsPointOnLS2D( poly_p2,p,inf_p)==1) 
		if( poly_p2(2,1) > poly_p1(2,1) ) 
		count=count+1; 
		end
	elseif (lsLSIntersect2D(p,inf_p,poly_p1,poly_p2)) 
		count=count+1; 
	end
end
if ( rem(count,2)==1 ) 
   in_flag=1;
   return;
else 
   in_flag=-1;
   return;
end