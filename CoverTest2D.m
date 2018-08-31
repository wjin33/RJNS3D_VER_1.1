function [rec_ij,lp_flag]=CoverTest2D(si,sj,sx,sy,line_p1,line_p2,colort)
%si,sj--the index of each square
%(sx,sy)--scale in 2 directions
%line_p1--3x1  starting point of the trace
%line_p2--3x1  ending point of the trace
%(x0,y0,z0)--coordinate of left lower point of covering points
%rec_ij--1x2, the index of each square
rec_ij=[];
x0=(si-1)*sx;  y0=(sj-1)*sy;
r_nvm=[x0, x0+sx, x0+sx, x0,    x0;
       y0, y0,    y0+sy, y0+sy, y0;
       0,  0,     0,     0,     0];
[row,lp_flag]=GetPolyLSIP2D(line_p1,line_p2,r_nvm);

if lp_flag~=-1 %draw the square
	x=[x0,x0+sx,x0+sx,x0+sx,x0+sx,x0,x0,x0];
	y=[y0,y0,y0,y0+sy,y0+sy,y0+sy,y0+sy,y0];
	line(x,y,'Color',colort,'LineWidth',2)
	rec_ij=[si,sj];
end
end