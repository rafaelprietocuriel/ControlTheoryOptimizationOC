function out=jumpdependentvar(hocTrj)

arcpos=arcposition(hocTrj);
out=zeros(size(hocTrj.y,1),2*jumpnum(hocTrj));
out(:,1:2)=hocTrj.y(:,1:2);
for ii=1:size(arcpos,2)-1
    out(:,2*ii+1:2*(ii+1))=hocTrj.y(:,arcpos(2,ii):arcpos(1,ii+1));
end
out(:,end-1:end)=hocTrj.y(:,end-1:end);