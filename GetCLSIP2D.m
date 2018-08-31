function [insc_p,is_flag]=GetCLSIP2D(cc,radius,poly_p1,poly_p2)
%poly_p1,poly_p2--3x1 matrix
%cc--3x1 matrix, the circle center 
%radius--radius of the joint
%is_flag--1:one point of intersection,  0:non-intersected,  2: two points of intersection
%insc_p--3x1 matrix or 3x2 matrix
insc_p=[];
is_flag=0;
u=(-(poly_p1(1,1)-cc(1,1))*(poly_p2(1,1)-poly_p1(1,1))-(poly_p1(2,1)-cc(2,1))*(poly_p2(2,1)-poly_p1(2,1)))/((poly_p2(1,1)-poly_p1(1,1))^2+(poly_p2(2,1)-poly_p1(2,1))^2);
dist=sqrt((poly_p1(1,1)+u*(poly_p2(1,1)-poly_p1(1,1))-cc(1,1))^2+(poly_p1(2,1)+u*(poly_p2(2,1)-poly_p1(2,1))-cc(2,1))^2);
if abs(dist-radius)<=1e-5
	if abs(u)<=1e-5
	u=0;
	elseif abs(u-1)<=1e-5
	u=1;
    end
    
	if (u>=0)&(u<=1)
		is_flag=1;
		insc_p=[poly_p1(1,1)+u*(poly_p2(1,1)-poly_p1(1,1));
                poly_p1(2,1)+u*(poly_p2(2,1)-poly_p1(2,1));
                0];	
		return
	end
elseif radius<dist
	is_flag=0;
	return
else
	u1=1/2/(poly_p2(1,1)^2-2*poly_p2(1,1)*poly_p1(1,1)+poly_p1(1,1)^2+poly_p2(2,1)^2-2*poly_p2(2,1)*poly_p1(2,1)+poly_p1(2,1)^2)*(2*cc(1,1)*poly_p2(1,1)-2*poly_p2(2,1)*poly_p1(2,1)+2*poly_p1(2,1)^2-2*poly_p2(1,1)*poly_p1(1,1)+2*poly_p1(1,1)^2-2*poly_p1(2,1)*cc(2,1)-2*poly_p1(1,1)*cc(1,1)+2*cc(2,1)*poly_p2(2,1)+2*(2*poly_p2(1,1)^2*poly_p1(2,1)*cc(2,1)-poly_p1(2,1)^2*cc(1,1)^2+poly_p1(2,1)^2*radius^2-poly_p1(1,1)^2*cc(2,1)^2+poly_p1(1,1)^2*radius^2-2*cc(1,1)*poly_p2(1,1)*poly_p1(2,1)*cc(2,1)-2*poly_p2(2,1)*poly_p1(2,1)*poly_p1(1,1)*cc(1,1)-2*poly_p2(1,1)*poly_p1(1,1)*poly_p1(2,1)*cc(2,1)+2*poly_p1(2,1)*cc(2,1)*poly_p1(1,1)*cc(1,1)-2*cc(2,1)*poly_p2(2,1)*poly_p1(1,1)*cc(1,1)+2*cc(1,1)*poly_p2(1,1)*poly_p1(2,1)^2+2*poly_p2(2,1)*poly_p1(2,1)*cc(1,1)^2-2*poly_p2(2,1)*poly_p1(2,1)*radius^2+2*poly_p2(1,1)*poly_p1(1,1)*cc(2,1)^2-2*poly_p2(1,1)*poly_p1(1,1)*radius^2+2*cc(2,1)*poly_p2(2,1)*poly_p1(1,1)^2-2*cc(1,1)*poly_p2(1,1)*poly_p2(2,1)*poly_p1(2,1)+2*cc(1,1)*poly_p2(1,1)*cc(2,1)*poly_p2(2,1)+2*poly_p2(2,1)*poly_p1(2,1)*poly_p2(1,1)*poly_p1(1,1)-2*poly_p2(1,1)*poly_p1(1,1)*cc(2,1)*poly_p2(2,1)-poly_p2(1,1)^2*poly_p1(2,1)^2+poly_p2(1,1)^2*radius^2-poly_p2(1,1)^2*cc(2,1)^2-poly_p2(2,1)^2*poly_p1(1,1)^2+poly_p2(2,1)^2*radius^2-poly_p2(2,1)^2*cc(1,1)^2+2*poly_p2(2,1)^2*poly_p1(1,1)*cc(1,1))^(1/2));
	u2=1/2/(poly_p2(1,1)^2-2*poly_p2(1,1)*poly_p1(1,1)+poly_p1(1,1)^2+poly_p2(2,1)^2-2*poly_p2(2,1)*poly_p1(2,1)+poly_p1(2,1)^2)*(2*cc(1,1)*poly_p2(1,1)-2*poly_p2(2,1)*poly_p1(2,1)+2*poly_p1(2,1)^2-2*poly_p2(1,1)*poly_p1(1,1)+2*poly_p1(1,1)^2-2*poly_p1(2,1)*cc(2,1)-2*poly_p1(1,1)*cc(1,1)+2*cc(2,1)*poly_p2(2,1)-2*(2*poly_p2(1,1)^2*poly_p1(2,1)*cc(2,1)-poly_p1(2,1)^2*cc(1,1)^2+poly_p1(2,1)^2*radius^2-poly_p1(1,1)^2*cc(2,1)^2+poly_p1(1,1)^2*radius^2-2*cc(1,1)*poly_p2(1,1)*poly_p1(2,1)*cc(2,1)-2*poly_p2(2,1)*poly_p1(2,1)*poly_p1(1,1)*cc(1,1)-2*poly_p2(1,1)*poly_p1(1,1)*poly_p1(2,1)*cc(2,1)+2*poly_p1(2,1)*cc(2,1)*poly_p1(1,1)*cc(1,1)-2*cc(2,1)*poly_p2(2,1)*poly_p1(1,1)*cc(1,1)+2*cc(1,1)*poly_p2(1,1)*poly_p1(2,1)^2+2*poly_p2(2,1)*poly_p1(2,1)*cc(1,1)^2-2*poly_p2(2,1)*poly_p1(2,1)*radius^2+2*poly_p2(1,1)*poly_p1(1,1)*cc(2,1)^2-2*poly_p2(1,1)*poly_p1(1,1)*radius^2+2*cc(2,1)*poly_p2(2,1)*poly_p1(1,1)^2-2*cc(1,1)*poly_p2(1,1)*poly_p2(2,1)*poly_p1(2,1)+2*cc(1,1)*poly_p2(1,1)*cc(2,1)*poly_p2(2,1)+2*poly_p2(2,1)*poly_p1(2,1)*poly_p2(1,1)*poly_p1(1,1)-2*poly_p2(1,1)*poly_p1(1,1)*cc(2,1)*poly_p2(2,1)-poly_p2(1,1)^2*poly_p1(2,1)^2+poly_p2(1,1)^2*radius^2-poly_p2(1,1)^2*cc(2,1)^2-poly_p2(2,1)^2*poly_p1(1,1)^2+poly_p2(2,1)^2*radius^2-poly_p2(2,1)^2*cc(1,1)^2+2*poly_p2(2,1)^2*poly_p1(1,1)*cc(1,1))^(1/2));
	count=0;
	if abs(u1)<=1e-5
	u1=0;
	elseif abs(u1-1)<=1e-5
	u1=1;
	end
	if abs(u2)<=1e-5
	u2=0;
	elseif abs(u2-1)<=1e-5
	u2=1;
	end	
	
	if (u1>=0)&(u1<=1)
		count=count+1;
	end
	if (u2>=0)&(u2<=1)
		count=count+1;
    end
    
	if count==0
		is_flag=0;
		return
	elseif count==1
		is_flag=1;
		if (u1>=0)&(u1<=1)
		insc_p=[poly_p1(1,1)+u1*(poly_p2(1,1)-poly_p1(1,1));
                poly_p1(2,1)+u1*(poly_p2(2,1)-poly_p1(2,1));
                0];
		end
		if (u2>=0)&(u2<=1)
		insc_p=[poly_p1(1,1)+u2*(poly_p2(1,1)-poly_p1(1,1));
                poly_p1(2,1)+u2*(poly_p2(2,1)-poly_p1(2,1));
                0];
		end
      	return
	else
		is_flag=2;
		insc_p=[poly_p1(1,1)+u1*(poly_p2(1,1)-poly_p1(1,1)),poly_p1(1,1)+u2*(poly_p2(1,1)-poly_p1(1,1));
                poly_p1(2,1)+u1*(poly_p2(2,1)-poly_p1(2,1)),poly_p1(2,1)+u2*(poly_p2(2,1)-poly_p1(2,1));
                0,0];
		return
	end
end