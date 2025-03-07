function ocgTrj=mergearc(ocgTrj,idx,varargin)

if isempty(ocgTrj) || isempty(idx)
    return
end
odenum=odenumber(ocgTrj);
maxodenum=max(odenum);
arcn=arcnum(ocgTrj);
arcpos=arcposition(ocgTrj);

ocTrj=struct(ocgTrj.octrajectory);
y=zeros(maxodenum,arcpos(2,arcn));
for ii=1:arcn
    y(1:odenum(ii),arcpos(1,ii):arcpos(2,ii))=ocTrj.y{ii};
end
ocTrj.y=y;
odenum(idx)=[];
ocgTrj=ocgtrajectory(mergearc(ocTrj,idx),odenum);