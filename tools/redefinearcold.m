function [sol,newarcargidx]=redefinearc(sol,newposition,arcid,opt)
%
% REDEFINEARC adds a new arc at a given position.
%
% REDEFINEARC(SOL,NEWPOSITION,ARCID) adds the new arc at NEWPOSITION with
% the arc identifier ARCID. The new arc has to be positioned in the
% interior
%
%
% [sol,newarcargidx]=redefinearcnew(sol,newposition,arcid,opt);
% return
newarcargidx=[];
if numel(newposition)==1
    newposition=[newposition;newposition];
end
newposition=newposition(:);
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];
numarc=numel(rightarcindex);
arcposition0=[leftarcindex;rightarcindex];

absx=transform2absx(sol);
if isfield(sol,'solverinfo') && isfield(sol.solverinfo,'tangent') && ~isempty(sol.solverinfo.tangent)
    tangent=sol.solverinfo.tangent;
    tangentend=tangent(numel(sol.y)+1:end);
    tangent=tangent(1:numel(sol.y));
    try
        tangent=reshape(tangent,size(sol.y));
    catch
        tangent=[];
    end
else
    tangent=[];
end
arcarg=zeros(1,numarc+1);
counter=0;
for ii=1:numarc
    newposcounter=find(newposition(1)>=leftarcindex(ii)&newposition(2)<=rightarcindex(ii),1);
    if ~isempty(newposcounter)
        if leftarcindex(ii)==newposition(1)
            if rightarcindex(ii)==newposition(2)
                counter=counter+1;
                arcarg(counter)=arcid;
            elseif rightarcindex(ii)>newposition(2)
                counter=counter+1;
                arcarg(counter)=arcid;
            end
        elseif leftarcindex(ii)<newposition(1) && newposition(2)<rightarcindex(ii)
            if rightarcindex(ii)==newposition(2)
                counter=counter+1;
                arcarg(counter)=sol.arcarg(ii);
                counter=counter+1;
                arcarg(counter)=arcid;
            elseif rightarcindex(ii)>newposition(2)
                counter=counter+1;
                arcarg(counter)=sol.arcarg(ii);
                counter=counter+1;
                arcarg(counter)=arcid;
            end
        elseif leftarcindex(ii)<newposition(1) && newposition(2)==rightarcindex(ii)
            counter=counter+1;
            arcarg(counter)=sol.arcarg(ii);
            counter=counter+1;
            arcarg(counter)=arcid;
        end
        if rightarcindex(ii)>newposition(2)
            counter=counter+1;
            arcarg(counter)=sol.arcarg(ii);
        end
    else
        counter=counter+1;
        arcarg(counter)=sol.arcarg(ii);
    end
end

if counter==numarc
    ocmatmsg('No new arc is added.\n')
    return
end
doublepos=unique(sort(newposition(:))).';
if ~(newposition(1,1)==1 && newposition(2,1)==1)
    doublepos(doublepos==1|doublepos==length(sol.x))=[];
    doublepos(ismember(doublepos,arcposition0(:)))=[];
else
    doublepos(doublepos==length(sol.x))=[];
end
if diff(newposition)==0 && ~(newposition(1,1)==1)
    doublepos=repmat(doublepos,1,2);
end
%doublepos((diff([-1 doublepos])==1))=[];
totalpos=sort([1:length(sol.x) doublepos(:).']);
absx=absx(totalpos);
sol.x=sol.x(totalpos);
sol.y=sol.y(:,totalpos);
arcposition=find(diff(sol.x)==0);
leftarcindex=[1 arcposition+1];
rightarcindex=[arcposition numel(sol.x)];
numarcold=numarc;
numarc=numel(rightarcindex);

sol.arcarg=arcarg;
sol.arcinterval=[absx(leftarcindex) absx(rightarcindex(end))];
idx=find(diff(leftarcindex)==1);
if length(idx)==1
    totalpos=sort([1:length(sol.x) leftarcindex([idx idx])]);
    sol.x=sol.x(totalpos);sol.y=sol.y(:,totalpos);
    rightarcindex(idx:end)=rightarcindex(idx:end)+2;leftarcindex(idx+1:end)=leftarcindex(idx+1:end)+2;
end
sol.arcposition=[leftarcindex;rightarcindex];
for ii=1:numarc
    sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii))=ii-1+transform2unit(sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii)));
end
if ~isempty(tangent)
    try
        tangent=tangent(:,totalpos);
        tangent=[tangent(:);zeros(numarc-numarcold,1);tangentend];
    catch
        tangent=[];
    end
end
sol.solverinfo.tangent=tangent(:);

if isfield(sol.solverinfo,'conttype')
    switch sol.solverinfo.conttype
        case 'extremal2ep'
        case 'extremale2ep'
            if isfield(solverinfo0,'jumpcostatecoord')
                sol.solverinfo.jumpold=solverinfo0.parameters(solverinfo0.jumpcostatecoord);
                sol.solverinfo.jumpcostatecoord=[];
            end
            sol.solverinfo.exactendpointcoord=solverinfo0.exactendpointcoord;
            sol.solverinfo.switchtimecoord=solverinfo0.exactendpointcoord(end)+[1:length(sol.arcinterval)-1];
            sol.solverinfo.parameters=[solverinfo0.parameters(solverinfo0.exactendpointcoord).' sol.arcinterval(2:end)];
        otherwise
            sol.solverinfo.parameters=[];
    end
else
    sol.solverinfo=[];
    %sol.solverinfo.parameters=[];
end
sol=mergearc(sol);

function sol=mergearc(sol)

idx=find(diff(sol.arcarg)==0);
if isempty(idx)
    removex=sol.arcposition(1,idx+1);
    sol.x(removex)=[];
    sol.y(:,removex)=[];
end
sol.arcarg(idx)=[];
sol.arcinterval(idx+1)=[];
arcposition=find(diff(sol.x)==0);
sol.arcposition=[1 arcposition+1;arcposition numel(sol.x)];
numarc=length(sol.arcposition(1,:));
for ii=1:numarc
    sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii))=ii-1+transform2unit(sol.x(sol.arcposition(1,ii):sol.arcposition(2,ii)));
end

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