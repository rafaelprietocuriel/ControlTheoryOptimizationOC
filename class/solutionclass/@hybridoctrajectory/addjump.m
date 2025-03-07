function hocTrj=addjump(hocTrj,jumptime,jumpid,varargin)

newposition=[];

if nargin>=4
    newposition=varargin{1};
end
indepvar=independentvar(hocTrj);
arcpos=arcposition(hocTrj);
arcint=arcinterval(hocTrj);
x0=initialtime(hocTrj);
arcn=arcnum(hocTrj);
diffarcinterval=diff(arcint);
t=x0;
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    t=[t t(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii)];
end
t=[t t(end)];

if isempty(newposition)
    if abs(jumptime-arcint(end))<1e-5
        hocTrj.jumparg(end)=jumpid;
        return
    elseif abs(jumptime-arcint(1))<1e-5
        hocTrj.jumparg(1)=jumpid;
        return
    else
        newposition=cont2idx(t,jumptime);
        if isempty(newposition)
            ocmatmsg('No jump is added.\n')
            return
        end
    end
end
newposition=newposition(1);
leftarcindex=arcpos(1,:);
rightarcindex=arcpos(2,:);
newposcounter=find(newposition>leftarcindex&newposition<rightarcindex,1);
if isempty(newposcounter)
    ocmatmsg('Jump occurs at jumptime.\n')
    return
end
sol=struct(hocTrj);
jumparg=zeros(1,length(sol.jumparg)+1);
jumparg([1 end])=sol.jumparg([1 end]);
arcarg=zeros(1,length(sol.arcarg)+1);
counter=0;
for ii=1:arcn
    if newposcounter==ii
        counter=counter+1;
        arcarg(counter)=sol.arcarg(ii);
        jumparg(counter)=sol.jumparg(ii);
        counter=counter+1;
        arcarg(counter)=sol.arcarg(ii);
        jumparg(counter)=jumpid;
    else
        counter=counter+1;
        arcarg(counter)=sol.arcarg(ii);
        jumparg(counter)=sol.jumparg(ii);
    end
end
totalpos=sort([1:length(sol.x) newposition]);
t=t(totalpos);
sol.x=sol.x(totalpos);
sol.y=sol.y(:,totalpos);
arcpos=find(diff(sol.x(2:end-1))==0)+1;
leftarcindex=[2 arcpos+1];
rightarcindex=[arcpos numel(sol.x)-1];
hocTrj.x=sol.x;
hocTrj.y=sol.y;
hocTrj.arcposition=[leftarcindex;rightarcindex];
arcn=length(hocTrj.arcposition(1,:));
for ii=1:arcn
    hocTrj.x(hocTrj.arcposition(1,ii):hocTrj.arcposition(2,ii))=ii-1+transform2unit(hocTrj.x(hocTrj.arcposition(1,ii):hocTrj.arcposition(2,ii)));
end
hocTrj.x([1 end])=hocTrj.x([2 end-1]);
hocTrj.arcarg=arcarg;
hocTrj.jumparg=jumparg;
hocTrj.arcinterval=[t(leftarcindex) t(end)];

%--------------------------------------------------------------------------
function xnew=transform2unit(x)
xnew=[];
n=numel(x);
if n<2
    return
end
l=x(n)-x(1);
if l==0
    xnew=linspace(0,1,n);
    return
end
xnew=(x-x(1))/l;