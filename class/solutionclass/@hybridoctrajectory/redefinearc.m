function ocTrj=redefinearc(ocTrj,newposition,arcid,opt)

if nargin==3
    opt=defaultocoptions;
end
structTrj=struct(ocTrj);
Xbd=structTrj.y(:,[1 end]);
structTrj.y(:,[1 end])=[];
structTrj.x([1 end])=[];
structTrj.arcposition=structTrj.arcposition-1;
newposition=newposition-1;

numarc=length(structTrj.arcarg);
jumparg=structTrj.jumparg;

[structTrj,newarcargidx]=redefinearc(structTrj,newposition,arcid,opt);
newjumparg=zeros(1,length(structTrj.arcarg)+1);
newjumparg(1)=jumparg(1);
for ii=1:numarc
    newjumparg(find(newarcargidx==ii,1,'last')+1)=jumparg(ii+1);
end
structTrj.jumparg=newjumparg;
structTrj.y=[Xbd(:,1) structTrj.y Xbd(:,2)];
structTrj.x=structTrj.x([1 1:end end]);
structTrj.arcposition=structTrj.arcposition+1;
ocTrj=hybridoctrajectory(structTrj);
