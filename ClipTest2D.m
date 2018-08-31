function [flag,u1,u2]=ClipTest2D(p,q,u1,u2)
%flag--variable identifier£¬ 0£ºabandoned£¬1£ºvisible
flag=1;
if(p<0.0)
	r=q/p;
	if (r>u2) 
		flag=0;
	elseif (r>u1)
		u1=r;     
      end
elseif(p>0.0)
	r=q/p;
	if (r<u1) 
		flag=0;
	elseif (r<u2)
		u2=r;
	end
elseif (q<0.0) 
	flag=0;
end