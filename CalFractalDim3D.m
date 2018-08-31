function rec_cv=CalFractalDim3D(sx0,sy0,sz0,mwh,mwv,ms,hm,cpm,occm,rm)
%rec_cv--1xn, the number of  box which cover the joints relat to the scale
%sx0,sy0,sz0--initial point of the model
%mwv--half length of the vertical axis of inner measure window
%mwh--half length of the horizontal axis of inner measure window
%ms--1xn measurement scale matrix 
%hm--1xn handle of each joint
%cpm--3xn center of each joint
%occm--2xn strike of each joint
%rm--1xn radius of each joint
%sqrcover--nxnxn used to record the box number

rec_cv=zeros(1,length(ms));
for m=1:length(ms)	
	hdiv=2*mwh/ms(m); vdiv=2*mwv/ms(m); wdiv=2*mwh/ms(m);
	sqrcover=zeros(hdiv,vdiv,wdiv);
	[row,col]=size(hm);
	for n=1:col
			vm=get(hm(n),'Vertices');
			vm=vm';
			temp=min(vm(1,:));
			%find the starting square in x-direction
			for rs_i=1:hdiv
				rs=sx0+(rs_i-1)*ms(m); re=sx0+rs_i*ms(m);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			%find the ending square in x-direction
			temp=max(vm(1,:));
			for re_i=rs_i:hdiv
				rs=sx0+(re_i-1)*ms(m); re=sx0+re_i*ms(m);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			%find the starting square in y-direction
			temp=min(vm(2,:));
			for rs_j=1:vdiv
				rs=sy0+(rs_j-1)*ms(m); re=sy0+rs_j*ms(m);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			%find the ending square in y-direction
			temp=max(vm(2,:));
			for re_j=rs_j:vdiv
				rs=sy0+(re_j-1)*ms(m); re=sy0+re_j*ms(m);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			%find the starting square in z-direction
			temp=min(vm(3,:));
			for rs_k=1:wdiv
				rs=sz0+(rs_k-1)*ms(m); re=sz0+rs_k*ms(m);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			%find the ending square in y-direction
			temp=max(vm(3,:));
			for re_k=rs_k:wdiv
				rs=sz0+(re_k-1)*ms(m); re=sz0+re_k*ms(m);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
            end	
            
            
			%draw the covering box
			hold on
            colort=rand(1,3);
			for i=rs_i:re_i
				for j=rs_j:re_j
					for k=rs_k:re_k	
						cv_flag=CoverTest3D(sx0,sy0,sz0,i,j,k,ms(m),ms(m),ms(m),hm(n),cpm(:,n),occm(:,n),rm(n));
						if cv_flag==1
							GenCuboidDis(sx0,sy0,sz0,i,j,k,ms(m),ms(m),ms(m),colort)
							sqrcover(i,j,k)=1;
						end
					end	
				end	
			end
	end
	
	covercounts=0;
	for n=1:wdiv
		temp=sqrcover(:,:,n);
		covercounts=covercounts+sum((sum(temp,1)));
	end
	rec_cv(m)=covercounts;
end
end