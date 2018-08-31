function es_lumtav=CalJointVolDens(jg_dd_dip,mw_dd_dip,es_lumtaa,m)
%jg_dd_dip--2x1 mean strike of the joint group, unit: radian
%mw_dd_dip--2x1 mean strike of the MW, unit: radian
%p1--3x1 mean normal vector of the joint group
%p2--3x1 normal vector of MW
%es_lumtav--volume density of the joint group 
%es_lumtaa--area density of the joint group 
%m--mean diameter of the joint disk 


alpha=jg_dd_dip(1,1); beta=jg_dd_dip(2,1);
%set the spherical coordinates parameter, theta and phi
theta=beta;
if alpha<=0.5*pi
	phi=0.5*pi-alpha;
else
	phi=2.5*pi-alpha;	
end
p1=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta);]; p1=PruneMartix(p1,1e-5);


alpha=mw_dd_dip(1,1); beta=mw_dd_dip(2,1);
%set the spherical coordinates parameter, theta and phi
theta=beta;
if alpha<=0.5*pi
	phi=0.5*pi-alpha;
else
	phi=2.5*pi-alpha;	
end
p2=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta);]; p2=PruneMartix(p2,1e-5);


%volume density=area density/(mean diameter*Sin(angle between MW and mean joint plane))
es_lumtav=es_lumtaa/(m*(1-dot(p1,p2)^2)^0.5);

end
