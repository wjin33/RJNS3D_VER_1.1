function [rec_cv,rec_cc]=CalFractalDim2D(mwh,mwv,ms,ipcm)

%rec_cv--1xn, the number of  square which cover the joints relat to the scale
%rec_cc--1xn, the number of joints covered relate to scale
%mwv--half length of the vertical axis of  measure window
%mwh--half length of the horizontal axis of  measure window
%ms--1xn measurement scale matrix
%ipcm--6xn the coordinates of each joint endpoints


rec_cv=zeros(1,length(ms));	rec_cc=zeros(1,length(ms));

for n=1:length(ms)	
	hdiv=2*mwh/ms(n); vdiv=2*mwv/ms(n); 
	sqrcover=zeros(hdiv,vdiv);%each cell equal 1 if the square cover joints
	sqrcount=zeros(hdiv,vdiv);%each cell store the joint number in this square
	ipcm=PruneMartix(ipcm,1e-5);
	[row,col]=size(ipcm);
	colort=rand(1,3);
	for k=1:col
			temp=min([ipcm(1,k),ipcm(4,k)]);
			%find the starting square in x-direction
			for rs_i=1:hdiv
				rs=(rs_i-1)*ms(n); 
                re=rs_i*ms(n);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
            end
			%find the ending square in x-direction
			temp=max([ipcm(1,k),ipcm(4,k)]);
			for re_i=rs_i:hdiv
				rs=(re_i-1)*ms(n); re=re_i*ms(n);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			
			%find the starting square in y-direction
			temp=min([ipcm(2,k),ipcm(5,k)]);
			for rs_j=1:vdiv
				rs=(rs_j-1)*ms(n); re=rs_j*ms(n);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			%find the starting square in y-direction
			temp=max([ipcm(2,k),ipcm(5,k)]);
			for re_j=rs_j:vdiv
				rs=(re_j-1)*ms(n); re=re_j*ms(n);
				if ((temp>rs)||(abs((temp-rs))<1e-5))&&((temp<re)||(abs((temp-re))<1e-5))
					break;
				end
			end
			
			%draw the covering square
			hold on
			for i=rs_i:re_i
				for j=rs_j:re_j
					[rec,lp_flag]=CoverTest2D(i,j,ms(n),ms(n),ipcm(1:3,k),ipcm(4:6,k),colort);
						if ~isempty(rec)
							sqrcover(rec(1),rec(2))=1;
							sqrcount(rec(1),rec(2))=sqrcount(rec(1),rec(2))+sqrcover(rec(1),rec(2));
						end
				end	
			end
	end
	
	rec_cv(1,n)=sum((sum(sqrcover,1)));
	rec_cc(1,n)=sum((sum(sqrcount,1)));

end

end
