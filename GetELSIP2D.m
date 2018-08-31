function [insc_p,is_flag]=GetELSIP2D(cc,major_axis,k_axis,poly_p1,poly_p2)
%poly_p1,poly_p2--2x1 matrix
%cc--2x1 matrix, the circle center 
%major_axis--the characteristic dimension
%k_axis--the retio of major to minor axis of a ellipse
%is_flag--1:one point of intersection,  0:non-intersected,  2: two points of intersection
%insc_p--2x1 matrix or 2x2 matrix
insc_p=[];
Q=[];
is_flag=0;
a=major_axis/2;
b=major_axis/(2*k_axis);
if abs(poly_p2(2,1)-poly_p1(2,1))<=1e-5
    if abs(abs(poly_p1(2,1))-b)<=1e-5
        Q=sort([poly_p1(1,1),poly_p2(1,1),0]);
        if Q(1,2)==0
            is_flag=1;
            insc_p=[0;poly_p1(2,1);0];
            return
        end    
    elseif abs(poly_p1(2,1))-b>1e-5
        is_flag=0;
        return
    else abs(poly_p1(2,1))-b<-1e-5
        x1=a*sqrt(1-(poly_p1(2,1)/b)^2);
        x2=-a*sqrt(1-(poly_p1(2,1)/b)^2);
        Q=sort([poly_p1(1,1),poly_p2(1,1)]);
        if Q(1,1)<=x1&&x1<=Q(1,2)
            is_flag=is_flag+1;
            p=[x1;poly_p1(2,1);0];
            insc_p=[insc_p,p];
        end
        if Q(1,1)<=x2&&x2<=Q(1,2)
            is_flag=is_flag+1;
            p=[x2;poly_p1(2,1);0];
            insc_p=[insc_p,p];
        end
        return
    end
elseif abs(poly_p2(1,1)-poly_p1(1,1))<=1e-5
    if abs(abs(poly_p1(1,1))-a)<=1e-5
        Q=sort([poly_p1(2,1),poly_p2(2,1),0]);
        if Q(1,2)==0
            is_flag=1;
            insc_p=[poly_p1(1,1);0;0];
            return
        end    
    elseif abs(poly_p1(1,1))-a>1e-5
        is_flag=0;
        return
    elseif abs(poly_p1(1,1))-a<-1e-5
        y1=b*sqrt(1-(poly_p1(1,1)/a)^2);
        y2=-b*sqrt(1-(poly_p1(1,1)/a)^2);
        Q=sort([poly_p1(2,1),poly_p2(2,1)]);
        if Q(1,1)<=y1&&y1<=Q(1,2)
            is_flag=is_flag+1;
            p=[poly_p1(1,1);y1;0];
            insc_p=[insc_p,p];
        end
        if Q(1,1)<=y2&&y2<=Q(1,2)
            is_flag=is_flag+1;
            p=[poly_p1(1,1);y2;0];
            insc_p=[insc_p,p];
        end
        return
    end
else
    k=(poly_p2(2,1)-poly_p1(2,1))/(poly_p2(1,1)-poly_p1(1,1));
    c=-k*poly_p1(1,1)+poly_p1(2,1);
    h0=sqrt((a^2*k^2+b^2)/(1+k^2));
    dist=abs(cc(2,1)-poly_p1(2,1)-k*(cc(1,1)-poly_p1(1,1)))/sqrt(k^2+1);
    if dist-h0>1e-5
        is_flag=0;
        return
    elseif abs(dist-h0)<=1e-5
        x1=(c/abs(c))*h0*sqrt(k^2+1)*(2*a^2*k)/(a^2*k^2+b^2);
        y1=(c/abs(c))*h0*sqrt(k^2+1)*(-2*a^2*k^2/(a^2*k^2+b^2)+1);
        lamuda=(x1-poly_p1(1,1))/(poly_p2(1,1)-poly_p1(1,1));
        if lamuda>=0
            insc_p=[x1;y1;0];
            is_flag=1;
             return
        end
    else
        x1=(-c*k*a^2+sqrt(k^2*c^2*a^4-(a^2*k^2+b^2)*(a^2*c^2-a^2*b^2)))/(a^2*k^2+b^2);
        x2=(-c*k*a^2-sqrt(k^2*c^2*a^4-(a^2*k^2+b^2)*(a^2*c^2-a^2*b^2)))/(a^2*k^2+b^2);
        y1=c+k*x1;
        y2=c+k*x2;
        lamuda=(x1-poly_p1(1,1))/(poly_p2(1,1)-x1);
        if lamuda>=0
            is_flag=is_flag+1;
            p=[x1;y1;0];
            insc_p=[insc_p,p];
        end
        lamuda=(x2-poly_p1(1,1))/(poly_p2(1,1)-x2);
        if lamuda>=0
            is_flag=is_flag+1;
            p=[x2;y2;0];
            insc_p=[insc_p,p];
        end
        return
    end
end

end