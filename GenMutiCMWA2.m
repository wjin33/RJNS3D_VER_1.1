function [hmw2,hmw3,len_tc_tra,len_td_tra,len_tcd_tra,len_tct_tra]=GenMutiCMWA2(x0,y0,z0,dx,dy,dz,joint_c,j_dd_dip,hm,cmwa,cmra,mw_dd_dip,ratio_h,ratio_v)
%(x0,y0,z0)--coordinate of the first vertex of cuboid
%(dx,dy,dz)--length,width,hight of the cuboid
%joint_c--3xn matrix, center of each joint
%j_dd_dip--2xn matrix, trend and dip of each joint
%hm--1xn matrix, handle of each joint
%cmwa--3xn centra of the MWs
%cmra-1xn radius of the inner circle
%mw_dd_dip--the strike of measure window
%ratio_h--the ratio of horizontal axis over radius, 
%ratio_v--the ratio of vertical  axis over radius,
%p_cn1--3x1 matrix, the global coordinates of right lower point of outer measuring window
%p_cn2--3x1 matrix, the global coordinates of left upper point of outer measuring window
%plane_p--6xn matrix, coordinates of the endpoint of trace in local measure window coordinate system
%p_cn1_1--the global coordinates of the left lower point of the first
%measuring window
%mwdya--half length of the vertical axis of outer measure window
%mwdxa--half length of the horizontal axis of outer measure window

mwdya=ratio_h*cmra;
mwdxa=ratio_v*cmra;

hmw2=figure;%generate a new figure 3-dimension
colort=rand(1,3);
%subsurface, generate vertexes from the first, anticlockwise
x=[x0,x0+dx,x0+dx,x0+dx,x0+dx,x0,   x0,   x0];
y=[y0,y0,   y0,   y0+dy,y0+dy,y0+dy,y0+dy,y0];
z=[z0,z0,   z0,   z0,   z0,   z0,   z0,   z0];
line(x,y,z,'Color',colort)
hold on

%topsurface, generate vertexes from the first, anticlockwise
x=[x0,x0+dx,x0+dx,x0+dx,x0+dx,x0,x0,x0];
y=[y0,y0,y0,y0+dy,y0+dy,y0+dy,y0+dy,y0];
z=[z0+dz,z0+dz,z0+dz,z0+dz,z0+dz,z0+dz,z0+dz,z0+dz];
line(x,y,z,'Color',colort)
hold on

%generate ridges from the first point, anticlockwise
x=[x0,x0];
y=[y0,y0];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

x=[x0+dx,x0+dx];
y=[y0,y0];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

x=[x0+dx,x0+dx];
y=[y0+dy,y0+dy];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
hold on

x=[x0,x0];
y=[y0+dy,y0+dy];
z=[z0,z0+dz];
line(x,y,z,'Color',colort)
DisAxis([0;0;0],3);
view(3);
axis equal;
hold on
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%generate multiple MW in the interested cuboid%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set c_ inside the inner MW;  d_ cutted one end  by inner MW, save the cutted trace
%cd_ cutted one end  by inner MW, save the whole trace; ct_ cross the inner MW
%len--means lenght; t--means the last row added the MW number
len_tc_tra=[]; len_td_tra=[]; len_tct_tra=[];  len_tcd_tra=[]; 
hmw3=figure;%generate a new figure 3-dimension


%calcualte the global coordinates of the left lower point of the first measuring window
cs_o=cmwa(:,1);
alpha=pi/2-mw_dd_dip(1,1);
beta=mw_dd_dip(2,1);
nto=[cos(alpha)*cos(beta),-sin(alpha),cos(alpha)*sin(beta);sin(alpha)*cos(beta),cos(alpha),sin(alpha)*sin(beta);-sin(beta),0,cos(beta);];
%left lower point of the first measure windom in local coordinate system
np_cn1=[mwdxa(1);-mwdya(1);0];
p_cn1_1=nto*np_cn1+cs_o;



[row,c_mw]=size(cmwa);
for k=1:c_mw
	set(0,'CurrentFigure',hmw2)
    %transform to current local measure windom coordinate system
	cs_o=cmwa(:,k);
	alpha=pi/2-mw_dd_dip(1,1);
	beta=mw_dd_dip(2,1);
	nto=[cos(alpha)*cos(beta),-sin(alpha),cos(alpha)*sin(beta);sin(alpha)*cos(beta),cos(alpha),sin(alpha)*sin(beta);-sin(beta),0,cos(beta);];
	%angular point of measure windom in local coordinate system
    np_cn1=[mwdxa(k);-mwdya(k);0];
	np_cn2=[-mwdxa(k);mwdya(k);0];
	%the four angular point in local coordinate system
	p_cn1=nto*np_cn1+cs_o;
	p_cn1=PruneMartix(p_cn1,1e-5);
	p_cn2=nto*np_cn2+cs_o;
	p_cn2=PruneMartix(p_cn2,1e-5);
    %transform to global coordinate system
	nrec_cn1=np_cn1;
	nrec_cn2=[np_cn1(1,1);np_cn2(2,1);0];
	nrec_cn3=np_cn2;
	nrec_cn4=[np_cn2(1,1);np_cn1(2,1);0];
    %transform to global coordinate system
	orec_cn1=nto*nrec_cn1+cs_o; 
	orec_cn2=nto*nrec_cn2+cs_o; 
	orec_cn3=nto*nrec_cn3+cs_o; 
	orec_cn4=nto*nrec_cn4+cs_o; 
	%draw the outer measure windom in global system
	line([orec_cn1(1,1),orec_cn2(1,1)],[orec_cn1(2,1),orec_cn2(2,1)],[orec_cn1(3,1),orec_cn2(3,1)]);
	line([orec_cn2(1,1),orec_cn3(1,1)],[orec_cn2(2,1),orec_cn3(2,1)],[orec_cn2(3,1),orec_cn3(3,1)]);
	line([orec_cn3(1,1),orec_cn4(1,1)],[orec_cn3(2,1),orec_cn4(2,1)],[orec_cn3(3,1),orec_cn4(3,1)]);
	line([orec_cn4(1,1),orec_cn1(1,1)],[orec_cn4(2,1),orec_cn1(2,1)],[orec_cn4(3,1),orec_cn1(3,1)]);

    
	%store the coordinates of the inner MW
	angle=0:pi/80:2*pi;
	[row,col]=size(angle);
	nc_nvm=zeros(3,col);
	oc_nvm=zeros(3,col);
	for i=1:col
		nc_nvm(:,i)=[cmra(k)*cos(angle(i));cmra(k)*sin(angle(i));0.0];
		oc_nvm(:,i)=nto*nc_nvm(:,i)+cs_o; 
    end
	oc_nvm=PruneMartix(oc_nvm,1e-5);
	%draw the inner circle measure windom in global system
	for i=1:(col-1)
		line([oc_nvm(1,i),oc_nvm(1,(i+1))],[oc_nvm(2,i),oc_nvm(2,(i+1))],[oc_nvm(3,i),oc_nvm(3,(i+1))]);
	end
    %%
	
   %considering the joint, from here; joint_c,j_dd_dip,hm
	[row,col]=size(hm);
	o_plane_p=[];
	for i=1:col
		[~,np1,np2,is_jmw_flag]=IsJointMWInsect(joint_c(:,i),j_dd_dip(:,i),hm(i),p_cn1,p_cn2,mw_dd_dip);
        if is_jmw_flag==1
        %attention: the original point of local coordinates system is the left lower point in function IsJointMWInsect
        	op1=nto*np1+p_cn1;
            op2=nto*np2+p_cn1;
            o_plane_p=[o_plane_p,[op1;op2];];
        end
    end
    %%
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%generate multiple MW in a new 2-D Window     %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	set(0,'CurrentFigure',hmw3)
   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[row,col]=size(o_plane_p);
    if col==0
        mes=strcat('There is no joint in ',int2str(k));
        disp(mes)
        continue 
    end
    % transform to local coordinates which use p_cn1_1 as it original point
 	cs_o=p_cn1_1;
	alpha=pi/2-mw_dd_dip(1,1);
	beta=mw_dd_dip(2,1);
	otn=[cos(alpha)*cos(beta),sin(alpha)*cos(beta),-sin(beta);-sin(alpha),cos(alpha),0;cos(alpha)*sin(beta),sin(alpha)*sin(beta),cos(beta);];
	[row,col]=size(o_plane_p);
	plane_p=zeros(6,col);
	for i=1:col
		op1=o_plane_p(1:3,i);
		op2=o_plane_p(4:6,i);
		np1=otn*(op1-cs_o);
		np2=otn*(op2-cs_o);
		plane_p(1:3,i)=np1;
		plane_p(4:6,i)=np2;
	end
	%rotate 90 degree around z-axis
	plane_p=[plane_p(2,:);-plane_p(1,:);plane_p(3,:);plane_p(5,:);-plane_p(4,:);plane_p(6,:)];
    [row,col]=size(plane_p);
    disp('shit')
	for i=1:col
		line([plane_p(1,i),plane_p(4,i)],[plane_p(2,i),plane_p(5,i)]);
	end
    %%

    %inner MW
	[row,col]=size(angle);
	for i=1:col
		nc_nvm(:,i)=otn*(oc_nvm(:,i)-cs_o); 
    end
	nc_nvm=[nc_nvm(2,:);-nc_nvm(1,:);nc_nvm(3,:)];
    for i=1:col-1
		line([nc_nvm(1,i),nc_nvm(1,(i+1))],[nc_nvm(2,i),nc_nvm(2,(i+1))]);
	end
    
     %outer MW
    nrec_cn1=otn*(orec_cn1-cs_o); 
	nrec_cn2=otn*(orec_cn2-cs_o); 
	nrec_cn3=otn*(orec_cn3-cs_o); 
	nrec_cn4=otn*(orec_cn4-cs_o); 
	%rotate 90 degree around z-axis
    nrec_cn1=[nrec_cn1(2,:);-nrec_cn1(1,:);nrec_cn1(3,:)];
    nrec_cn2=[nrec_cn2(2,:);-nrec_cn2(1,:);nrec_cn2(3,:)];
    nrec_cn3=[nrec_cn3(2,:);-nrec_cn3(1,:);nrec_cn3(3,:)];
    nrec_cn4=[nrec_cn4(2,:);-nrec_cn4(1,:);nrec_cn4(3,:)];
	%censor
	plane_p=PruneMartix(plane_p,1e-5);
	nc_nvm=PruneMartix(nc_nvm,1e-5);
	nrec_cn1=PruneMartix(nrec_cn1,1e-5);
	nrec_cn2=PruneMartix(nrec_cn2,1e-5);
	nrec_cn3=PruneMartix(nrec_cn3,1e-5);
	nrec_cn4=PruneMartix(nrec_cn4,1e-5);
    %draw the outer measure windom in local system
	line([nrec_cn1(1,1),nrec_cn2(1,1)],[nrec_cn1(2,1),nrec_cn2(2,1)]);
	line([nrec_cn2(1,1),nrec_cn3(1,1)],[nrec_cn2(2,1),nrec_cn3(2,1)]);
	line([nrec_cn3(1,1),nrec_cn4(1,1)],[nrec_cn3(2,1),nrec_cn4(2,1)]);
	line([nrec_cn4(1,1),nrec_cn1(1,1)],[nrec_cn4(2,1),nrec_cn1(2,1)]);

    %%
    
	%set c_ inside the inner MW;  d_ cutted one end  by inner MW, save the cutted trace
    %cd_ cutted one end  by inner MW, save the whole trace; ct_ cross the inner MW
	c_tra=[]; d_tra=[]; cd_tra=[]; ct_tra=[];  
	len_c_tra=[]; len_d_tra=[]; len_cd_tra=[]; len_ct_tra=[];
	%c_tra; d_tra; cd_tra; ct_tra--temporary variable
	%estimate the intersection bewteen trace and inner(outer) MW
	[row,col]=size(plane_p);
	for i=1:col
        [rec_p,lp_flag]=GetPolyLSIP2D(plane_p(1:3,i),plane_p(4:6,i),nc_nvm);
		if lp_flag==0
			c_tra=[c_tra,plane_p(:,i)];
			temp=plane_p(:,i);
			lens=sqrt(((temp(1,1)-temp(4,1))^2+(temp(2,1)-temp(5,1))^2));
			len_c_tra=[len_c_tra,lens];
					
		elseif lp_flag==1
				d_tra=[d_tra,rec_p];
				temp=rec_p;
				lens=sqrt(((temp(1,1)-temp(4,1))^2+(temp(2,1)-temp(5,1))^2));
				len_d_tra=[len_d_tra,lens];
				
				cd_tra=[cd_tra,plane_p(:,i)];
				temp=plane_p(:,i);		
				lens=sqrt(((temp(1,1)-temp(4,1))^2+(temp(2,1)-temp(5,1))^2));
				len_cd_tra=[len_cd_tra,lens];
				%%%%

		elseif lp_flag==2
			ct_tra=[ct_tra,plane_p(:,i)];
			temp=plane_p(:,i);
			lens=sqrt(((temp(1,1)-temp(4,1))^2+(temp(2,1)-temp(5,1))^2));
			len_ct_tra=[len_ct_tra,lens];				
		else
			continue	
		end
	end

	
	[row,col]=size(len_c_tra);
	if col~=0
		c_tra=[c_tra;k*ones(1,col)];
		len_c_tra=[len_c_tra;k*ones(1,col)];
		len_tc_tra=[len_tc_tra,len_c_tra];
	else
		len_c_tra=[0;k];
		len_tc_tra=[len_tc_tra,len_c_tra];		
	end
		
	[row,col]=size(len_d_tra);
	if col~=0
		d_tra=[d_tra;k*ones(1,col)];
		len_d_tra=[len_d_tra;k*ones(1,col)];
		len_td_tra=[len_td_tra,len_d_tra];
	else
		len_d_tra=[0;k];
		len_td_tra=[len_td_tra,len_d_tra];			
	end
		
	[row,col]=size(len_cd_tra);
	if col~=0
		cd_tra=[cd_tra;k*ones(1,col)];
		len_cd_tra=[len_cd_tra;k*ones(1,col)];
		len_tcd_tra=[len_tcd_tra,len_cd_tra];	
	else
		len_cd_tra=[0;k];
		len_tcd_tra=[len_tcd_tra,len_cd_tra];			
	end
			
	[row,col]=size(len_ct_tra);
	if col~=0	
		ct_tra=[ct_tra;k*ones(1,col)];
		len_ct_tra=[len_ct_tra;k*ones(1,col)];
		len_tct_tra=[len_tct_tra,len_ct_tra];
	else
		len_ct_tra=[0;k];
		len_tct_tra=[len_tct_tra,len_ct_tra];			
    end

    [row,col]=size(ct_tra);
	if col~=0
		colort=[1,0,0];  %red
		for i=1:col
			line(([ct_tra(1,i),ct_tra(4,i)]),([ct_tra(2,i),ct_tra(5,i)]),'Color',colort,'LineWidth',1);
        end
		theta=atan(abs((ct_tra(5,1)-ct_tra(2,1))/(ct_tra(4,1)-ct_tra(1,1))));
	end
		
	[row,col]=size(cd_tra);
	if col~=0
		colort=[0,1,0];  %green
		for i=1:col
			line(([cd_tra(1,i),cd_tra(4,i)]),([cd_tra(2,i),cd_tra(5,i)]),'Color',colort,'LineWidth',1);
        end
		theta=atan(abs((cd_tra(5,1)-cd_tra(2,1))/(cd_tra(4,1)-cd_tra(1,1))));
	end
		
	[row,col]=size(c_tra);
	if col~=0
		colort=[0,1,1];  %cyan
		for i=1:col
			line(([c_tra(1,i),c_tra(4,i)]),([c_tra(2,i),c_tra(5,i)]),([c_tra(3,i),c_tra(6,i)]),'Color',colort,'LineWidth',1);
        end
		theta=atan(abs((c_tra(5,1)-c_tra(2,1))/(c_tra(4,1)-c_tra(1,1))));
    end
		
end

axis equal

end