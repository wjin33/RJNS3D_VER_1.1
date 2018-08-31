function DisAxis(pcen,axislen)
%pcen--3x1 matrix, coordinate of the original point
%axislen--scalar,length of the axis
hold on

%local axis-z, red
dv=[0.0,0.0,axislen]';
dv=dv+pcen;
xvd=[dv(1,1),pcen(1,1)];
yvd=[dv(2,1),pcen(2,1)];
zvd=[dv(3,1),pcen(3,1)];
line(xvd,yvd,zvd,'Color',[1,0,0],'LineWidth',4);

%local axis-x, green
dv=[axislen,0.0,0.0]';
dv=dv+pcen;
xvd=[dv(1,1),pcen(1,1)];
yvd=[dv(2,1),pcen(2,1)];
zvd=[dv(3,1),pcen(3,1)];
line(xvd,yvd,zvd,'Color',[0,1,0],'LineWidth',4);

%local axis-y, blue
dv=[0.0,axislen,0.0]';
dv=dv+pcen;
xvd=[dv(1,1),pcen(1,1)];
yvd=[dv(2,1),pcen(2,1)];
zvd=[dv(3,1),pcen(3,1)];
line(xvd,yvd,zvd,'Color',[0,0,1],'LineWidth',4);
hold on