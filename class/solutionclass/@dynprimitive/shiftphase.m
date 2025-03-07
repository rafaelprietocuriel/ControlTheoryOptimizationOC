function dynPrim=shiftphase(dynPrim,phaseidx)

if ~isperiodic(dynPrim) || phaseidx==1
    return
end
arcarg=arcargument(dynPrim);
arcpos=arcposition(dynPrim);
depvar=dependentvar(dynPrim);
indepvar=independentvar(dynPrim);
arcint=arcinterval(dynPrim);
if phaseidx==arcpos(2,end)
    return
end
diffarcinterval=diff(arcint);
s=initialtime(dynPrim);
numarc=length(arcarg);
for ii=1:numarc
    arcp=arcpos(1,ii):arcpos(2,ii);
    s=s(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii);
    x(1,arcp)=s;
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
    arcshift=find(phaseidx<arcpos(2,:),1,'first');
    dynPrimStruct.octrajectory.arcarg=arcarg([arcshift:numarc 2:arcshift]);

    newx=[x(newpos1)-x(newpos1(1)) x(newpos2)+x(newpos1(end))-x(newpos1(1))];
    newindepvar=[indepvar(newpos1)-indepvar(newpos1(1)) indepvar(newpos2)+indepvar(newpos1(end))-indepvar(newpos1(1))];

    l1=[arcpos(1,arcshift+1:end)-phaseidx+1 arcpos(1,2:arcshift)+arcpos(2,end)-phaseidx];
    newarcposition=[1 l1;l1-1 arcpos(2,end)];
    %idx=find(~diff(newindepvar));
    %newarcposition=[1 idx+1;idx length(newindepvar)];

    removearc=[];
    timeinterval=zeros(1,numarc);
    for ii=1:numarc
        timeinterval(ii)=newx(newarcposition(2,ii))-newx(newarcposition(1,ii));
        if timeinterval(ii)>0
        newindepvar(newarcposition(1,ii):newarcposition(2,ii))=ii-1+(newindepvar(newarcposition(1,ii):newarcposition(2,ii))-newindepvar(newarcposition(1,ii)))/(newindepvar(newarcposition(2,ii))-newindepvar(newarcposition(1,ii)));
        else
            removearc=[removearc ii];
        end
    end
    timeinterval(removearc)=[];
    newarcposition(:,removearc)=[];
    dynPrimStruct.octrajectory.arcarg(removearc)=[];
    dynPrimStruct.octrajectory.arcinterval=[initialtime(dynPrim) initialtime(dynPrim)+cumsum(timeinterval)];
    dynPrimStruct.octrajectory.x=newindepvar;
    dynPrimStruct.octrajectory.arcposition=newarcposition;
else
end
freepar=parameters(dynPrim);
if isfield(dynPrimStruct.octrajectory.solverinfo,'jumpcostatecoord')
    jump=freepar(dynPrimStruct.octrajectory.solverinfo.jumpcostatecoord);
else
    jump=[];
end
dynPrimStruct.octrajectory.solverinfo.parameters=[jump cumsum(timeinterval) freepar(end)];
dynPrimStruct.octrajectory.solverinfo.tangent=[];
dynPrimStruct.octrajectory.solverinfo.coeff=[];
dynPrimStruct.octrajectory.solverinfo.tmesh=[];
dynPrim=dynprimitive(dynPrimStruct);