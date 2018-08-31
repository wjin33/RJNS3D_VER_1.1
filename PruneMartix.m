function nm=PruneMartix(im,toler)
%im--input martix 
%toler--scalar
[row,col]=find(abs(im)<=toler);
[ri,ci]=size(row);
for i=1:ri
im(row(i,1),col(i,1))=0;
end
nm=im;