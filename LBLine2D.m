function [p1,p2,is_flag]=LBLine2D(xwmin,ywmin,xwmax,ywmax,x1,y1,x2,y2)
%is_flag--1:intersect£¬0:non-intersect
%p1,p2--3x1 matrix£¬the coordinates of  intersection point
u1=0.0; u2=1.0; %0£ºabandoned£¬1£ºvisible
p1=zeros(3,1);   
p2=zeros(3,1);
is_flag=0;
dx=x2-x1;
[flag,u1,u2]=ClipTest2D(-dx,x1-xwmin,u1,u2);
if (flag)
	[flag,u1,u2]=ClipTest2D(dx,xwmax-x1,u1,u2);
	if (flag)
		dy=y2-y1;
		[flag,u1,u2]=ClipTest2D(-dy,y1-ywmin,u1,u2);
		if (flag)
			[flag,u1,u2]=ClipTest2D(dy,ywmax-y1,u1,u2);
			if (flag)
				if (u2<1.0)
   				    	x2=x1+u2*dx;
					y2=y1+u2*dy;
				end
				if (u1>0.0)
					x1=x1+u1*dx;
					y1=y1+u1*dy;
				end
				is_flag=1;
				p1=[x1;y1;0];
				p2=[x2;y2;0];
			end    
		end
	end
end	