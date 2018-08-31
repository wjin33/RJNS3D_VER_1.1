function jt_dd_dip=TransJointDDDips(mw_strike,j_dd_dip)
%mw_strike--1x1 measure window trend, unit: radian
%mw_dip--dip of measure window: 90 degree
%j_dd_dip--2xn the original trend and dip of joint group
%jt_dd_dip--2xn the transformed trend and dip of joint group


alpha=pi/2-mw_strike; mw_dip=0.5*pi; beta=-(0.5*pi-mw_dip);
nto=[cos(alpha)*cos(beta),-sin(alpha),cos(alpha)*sin(beta);
     sin(alpha)*cos(beta),cos(alpha), sin(alpha)*sin(beta);
     -sin(beta),          0,          cos(beta);           ];
nto=PruneMartix(nto,1e-5);
%%transform to N-Y,E-X globle coordinates,j_dd are trend, j_dips are dip
[row,col]=size(j_dd_dip);
jt_dd_dip=zeros(row,col);
for j=1:col
	j_dd=j_dd_dip(1,j); j_dip=j_dd_dip(2,j);
	theta=j_dip;
	if j_dd<=0.5*pi
		phi=0.5*pi-j_dd;
	else
		phi=2.5*pi-j_dd;	
	end
	j_ov=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta);];
	j_nv=[(j_ov')*nto(:,1);(j_ov')*nto(:,2);(j_ov')*nto(:,3);];
	j_nv=PruneMartix(j_nv,1e-5);
	theta_t=acos(j_nv(3,1));
    dip_t=theta_t;
    
    
	cos_phi_t=j_nv(1,1)/sin(theta_t); sin_phi_t=j_nv(2,1)/sin(theta_t);
	if sin_phi_t>=0
		if cos_phi_t>=0
			phi_t=asin(sin_phi_t);	
		else
			phi_t=acos(cos_phi_t);
		end
	else
		if cos_phi_t<=0
			phi_t=pi+asin(-sin_phi_t);	
		else
			phi_t=2*pi+asin(sin_phi_t);
		end
    end
	if phi_t>0.5*pi
		dd_t=2.5*pi-phi_t;
	else
		dd_t=0.5*pi-phi_t;
    end
	jt_dd_dip(1,j)=dd_t;
	jt_dd_dip(2,j)=dip_t;

end

end
