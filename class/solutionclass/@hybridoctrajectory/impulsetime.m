function out=impulsetime(hocTrj)

jumparg=jumpargument(hocTrj);
jumparg([1 end])=-1;
out=arcinterval(hocTrj);

out(jumparg==0)=[];