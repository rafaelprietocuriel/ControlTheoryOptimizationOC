function dynPrim=shiftphasenew(dynPrim,phaseidx)

if ~isperiodic(dynPrim) || phaseidx==1
    return
end
arcarg=arcargument(dynPrim);
arcpos=arcposition(dynPrim);
depvar=dependentvar(dynPrim);
indepvar=independentvar(dynPrim);
arcint=arcinterval(dynPrim);
if numel(phaseidx)>1
    return
end
if phaseidx==arcpos(2,end) || phaseidx==1
    return
end
if any(arcpos(1,:)-phaseidx==0) || any(arcpos(2,:)-phaseidx==0)
    return
end

diffarcinterval=diff(arcint);
s=initialtime(dynPrim);
numarc=length(arcarg);
for ii=1:numarc
    arcp=arcpos(1,ii):arcpos(2,ii);
    s=s(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii);
    t(1,arcp)=s;
end
dynPrimStruct.octrajectory=struct(octrajectory(dynPrim));
dynPrimStruct.octrajectory.linearization=[];
dynPrimStruct.octrajectory.solver=[];
%dynPrimStruct.octrajectory.solverinfo=[];
dynPrimStruct.period=period(dynPrim);
if arcarg(1)==arcarg(end)
    newpos1=phaseidx:arcpos(2,end);
    newpos2=2:phaseidx;
    newpos=[newpos1 newpos2];
    dynPrimStruct.octrajectory.y=depvar(:,newpos);
    k=find(arcpos(1,:)-phaseidx<0,1,'last');
    narcarg=[arcarg(k:end-1) arcarg(1:k) ];
    tnew=[t(phaseidx:end) t(end)+t(2:phaseidx)];
    tnew=tnew-tnew(1);

    narcpos=find(diff(tnew)==0);
    narcpos=[1 narcpos+1;narcpos length(tnew)];

    for ii=1:numarc
        arcp=narcpos(1,ii):narcpos(2,ii);
        tarc=tnew(arcp);
        timeinterval(ii)=tarc(end)-tarc(1);
        s(1,arcp)=ii-1+(tarc-tarc(1))/timeinterval(ii);
    end
    %timeinterval=cumsum([0 timeinterval]);
    dynPrimStruct.octrajectory.arcarg=narcarg;
    dynPrimStruct.octrajectory.arcinterval=[initialtime(dynPrim) initialtime(dynPrim)+cumsum(timeinterval)];
    dynPrimStruct.octrajectory.x=s;
    dynPrimStruct.octrajectory.arcposition=narcpos;
else
end
%freepar=parameters(dynPrim);

%dynPrimStruct.octrajectory.solverinfo.parameters=[timeinterval(2:end) freepar(end)];
dynPrimStruct.octrajectory.solverinfo.parameters=[timeinterval(2:end) period(dynPrim)];
dynPrimStruct.octrajectory.solverinfo.tangent=[];
dynPrimStruct.octrajectory.solverinfo.coeff=[];
dynPrimStruct.octrajectory.solverinfo.tmesh=[];
dynPrim=dynprimitive(dynPrimStruct);