function [sol,newarcargindex]=redefinearc(sol,newposition,newarcid,opt)
%
% REDEFINEARC adds a new arc at a given position.
%
% REDEFINEARC(SOL,NEWPOSITION,ARCID) adds the new arc at NEWPOSITION with
% the arc identifier ARCID. The new arc has to be positioned in the
% interior
%
% Each column of NEWPOSITION corresponds to a new arc. If two or more new
% arcs are eqach new arc has to be contained in a different arc of the old
% arcs.
%
% If NEWPOSITION is a single value the corresponding arc is trivially
% discretized by the number of 'TrivialArcMeshNum' defined in the GENERAL
% options.
%
% REDEFINEARC(SOL,NEWPOSITION,ARCID,OPT) with OPT a different number of
% 'TrivialArcMeshNum' can be provided.
%
% [SOL,NEWARCARGINDEX]=REDEFINEARC(SOL,...) SOL contains the new solution
% structure. NEWARCARGINDEX returns the pattern, where 0 indicates the
% place, where an old arc identifier is located within the new arc
% identifiers.

newarcargindex=[];
if numel(newposition)==1
    newposition=[newposition;newposition];
end
if any(newposition(1,2:end)-newposition(2,1:end-1)<=0)
    ocmatmsg('Intervals of new positions must not intersect.\n')
    return
end
if size(newposition,2)-length(newarcid)~=0
    ocmatmsg('Number of new intervals and number of identifiers have to be same.\n')
    return
end
if any(diff(newposition)<0)
    ocmatmsg('Intervals of new positions have to be increasing.\n')
    return
end
if nargin==3
    opt=defaultocoptions;
end
initnewarcmeshnum=getocoptions(opt,'GENERAL','TrivialArcMeshNum');

arcposition0=sol.arcposition;
arcarg0=sol.arcarg;
arcinterval0=sol.arcinterval;
y0=sol.y;
x0=sol.x;
if isfield(sol.solverinfo,'meshvalcoord') && isfield(sol.solverinfo,'parametercoord') && isfield(sol.solverinfo,'switchtimecoord') && isfield(sol.solverinfo,'tangent') && ~isempty(sol.solverinfo.tangent)
    tgt0=sol.solverinfo.tangent;
    coeff0=sol.solverinfo.coeff;
    tgtpar0=tgt0(sol.solverinfo.parametercoord);
    coeffpar0=coeff0(sol.solverinfo.parametercoord);
    tgt0(sol.solverinfo.parametercoord)=[];
    tgt0=tgt0(sol.solverinfo.meshvalcoord);
    coeff0(sol.solverinfo.parametercoord)=[];
    coeff0=coeff0(sol.solverinfo.meshvalcoord);
    switchtimecoord0=sol.solverinfo.switchtimecoord;
else
    tgt0=[];
    tgtpar0=[];
    coeff0=[];
    coeffpar0=[];
    switchtimecoord0=[];
end

numarc=length(arcarg0);
if newposition(1,1)<1 || newposition(2,end)>arcposition0(2,end)
    ocmatmsg('New positions have to be positive and must not exceed positions of old interval.\n')
    return
end

absx=transform2absx(sol);
newarcintervaltime=reshape(absx(newposition),2,[]);
newarcarg=arcarg0;
newarcargindex=1:length(arcarg0);
arcargindex0=newarcargindex;
arcposition=arcposition0;
arcinterval=arcinterval0;
y=y0;
x=x0;
tgt=tgt0;
tgtpar=tgtpar0(:).';
coeff=coeff0;
coeffpar=coeffpar0(:).';
switchtimecoord=switchtimecoord0;
oldarcpos=[];
for ii=size(newposition,2):-1:1
    % assure that new intervals are entirely included in one of the
    % original intervals and
    % bdflag = -1 new interval starts at left side of an old interval
    % bdflag = 1 new interval starts at right side of an old interval
    % bdflag = 0 new interval starts in the interior of an old interval
    newarcpos=find(newposition(1,ii)-arcposition0(1,1:end)>=0,1,'last');
    if ~isempty(oldarcpos)
        if oldarcpos==newarcpos
            newarcargindex=[];
            ocmatmsg('New intervals must not be contained in the same old interval.\n')
            return
        end
    end
    oldarcpos=newarcpos;
    if newposition(2,ii)-arcposition0(2,newarcpos)>0
        ocmatmsg('Intervals of new positions have to be fully part of old arc.\n')
        return
    end
    if arcposition0(1,newarcpos)-newposition(1,ii)==0
        if newposition(2,ii)-arcposition0(2,newarcpos)==0
            ocmatmsg('Interval of new positions replaces an old arc. Use ''removearc'' or ''mergearc'' in this case.\n')
            return
        end
        bdflag=-1;
    elseif newposition(2,ii)-arcposition0(2,newarcpos)==0
        bdflag=1;
    else
        bdflag=0;
    end
    if newposition(2,ii)-newposition(1,ii)<=1
        newarcmeshnum=initnewarcmeshnum;
        if newposition(2,ii)==newposition(1,ii)
            ynew=repmat(y0(:,newposition(2,ii)),1,newarcmeshnum);
        else
            ynew=repmat(y0(:,newposition(1,ii)),1,newarcmeshnum)+repmat(linspace(0,1,newarcmeshnum),size(y0,1),1).*(repmat(y0(:,newposition(2,ii))-y0(:,newposition(1,ii)),1,newarcmeshnum));
        end
        if ~isempty(tgt)
            tgtnew=repmat(zeros(size(y0(:,newposition(2,ii)))),1,newarcmeshnum);
        end
    else
        newarcmeshnum=1;
        ynew=y0(:,newposition(1,ii):newposition(2,ii));
        if ~isempty(tgt)
            tgtnew=zeros(size(y0(:,newposition(1,ii):newposition(2,ii))));
        end
    end
    switch bdflag
        case 0
            newarcarg=[arcarg0(1:newarcpos) newarcid(ii) newarcarg(newarcpos:end)]; % new arcargument
            if ~isempty(tgt)
                tgtspar=[tgtpar(switchtimecoord(1:newarcpos-1)) -1 1 tgtpar(switchtimecoord(newarcpos+1:end))];
                coeffspar=[coeffpar(switchtimecoord(1:newarcpos-1)) x0([newposition(1,ii) newposition(2,ii)]) coeffpar(switchtimecoord(newarcpos+1:end))];
                if ~isempty(switchtimecoord)
                    tgtpar=[tgtpar(1:switchtimecoord(1)-1) tgtspar tgtpar(switchtimecoord(end)+1:end)];
                    coeffpar=[coeffpar(1:switchtimecoord(1)-1) coeffspar coeffpar(switchtimecoord(end)+1:end)];
                else
                    tgtpar=[tgtpar(1:end-1) tgtspar];
                    coeffpar=[coeffpar(1:end-1) coeffspar];
                end
                tgtpar=[tgtpar 0];
                coeffpar=[coeffpar coeffpar0(end)];
            end
            newarcargindex=[arcargindex0(1:newarcpos) NaN newarcargindex(newarcpos:end)]; % pattern for the "old" arc identifiers within the new arc identifiers
            if newarcmeshnum>1
                if newposition(2,ii)==newposition(1,ii) % a new arc is added at a single point
                    arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);newposition(1,ii)] [newposition(1,ii)+1;newarcmeshnum+1+newposition(1,ii)-1] newarcmeshnum+1+[newposition(1,ii);arcposition(2,newarcpos)] newarcmeshnum+1+arcposition(:,newarcpos+1:end)];
                    xnewl=x0(arcposition0(1,newarcpos):newposition(1,ii));
                    y=[y0(:,1:newposition(1,ii)) ynew y(:,newposition(2,ii):end)];
                else % a new arc is added between two adjacent points
                    arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);newposition(1,ii)] [newposition(1,ii)+1;newarcmeshnum+newposition(1,ii)] newarcmeshnum+[newposition(1,ii)+1;arcposition(2,newarcpos)] newarcmeshnum+arcposition(:,newarcpos+1:end)];
                    xnewl=x0(arcposition0(1,newarcpos):newposition(1,ii));
                    y=[y0(:,1:newposition(1,ii)) ynew y(:,newposition(2,ii):end)];
                end
                xnewm=linspace(0,1,newarcmeshnum)+newarcpos;
                xnewr=x0(newposition(2,ii):arcposition0(2,newarcpos));
            else % a new arc is added between more than two consecutive points
                arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);newposition(1,ii)] [newposition(1,ii)+1;newposition(2,ii)+1] [newposition(2,ii)+2;arcposition(2,newarcpos)+2] 2+arcposition(:,newarcpos+1:end)];
                xnewl=x0(arcposition0(1,newarcpos):newposition(1,ii));
                xnewm=x0(newposition(1,ii):newposition(2,ii));
                xnewr=x0(newposition(2,ii):arcposition0(2,newarcpos));
                y=[y0(:,1:newposition(1,ii)) ynew y(:,newposition(2,ii):end)];
            end
            if ~isempty(tgt)
                tgt=[tgt0(:,1:newposition(1,ii)) tgtnew tgt(:,newposition(2,ii):end)];
                coeff=y;
            end
            arcinterval=[arcinterval0(1:newarcpos) newarcintervaltime(:,ii).' arcinterval(newarcpos+1:end)];

            % normalize time intervals
            xnewl=(xnewl-xnewl(1))/(xnewl(end)-xnewl(1))-1+newarcpos;
            xnewm=(xnewm-xnewm(1))/(xnewm(end)-xnewm(1))+newarcpos;
            xnewr=(xnewr-xnewr(1))/(xnewr(end)-xnewr(1))+1+newarcpos;

            if newarcpos<numarc
                if newarcpos>1
                    x=[x0(:,1:arcposition0(1,newarcpos)-1) xnewl xnewm xnewr 2+x(:,arcposition0(1,newarcpos+1):end)];
                else
                    x=[xnewl xnewm xnewr 2+x(:,arcposition0(1,newarcpos+1):end)];
                end
            else
                x=[x0(:,1:arcposition0(1,newarcpos)-1) xnewl xnewm xnewr];
            end
            if ~isempty(tgt)
                parametercoord=numel(y)+(1:length(tgtpar));
                meshvalcoord=reshape(1:numel(y),size(y,1),[]);
                tgt=[tgt(:).' tgtpar];
                tgt=tgt/norm(tgt);
                coeff=[coeff(:);coeffpar(:)];
            end

        case -1
            newarcarg=[arcarg0(1:newarcpos-1) newarcid(ii) newarcarg(newarcpos:end)];
            newarcargindex=[arcargindex0(1:newarcpos-1) NaN newarcargindex(newarcpos:end)];
            xnewl=[];
            if newarcmeshnum>1
                xnewm=linspace(0,1,newarcmeshnum);
                if newposition(2,ii)==newposition(1,ii)
                    arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);arcposition0(1,newarcpos)+newarcmeshnum-1] newarcmeshnum+[newposition(2,ii);arcposition(2,newarcpos)] newarcmeshnum+arcposition(:,newarcpos+1:end)];
                else
                    arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);arcposition0(1,newarcpos)+newarcmeshnum-1] newarcmeshnum+[newposition(2,ii)-1;arcposition(2,newarcpos)-1] newarcmeshnum-1+arcposition(:,newarcpos+1:end)];
                end
            else
                xnewm=x0(newposition(1,ii):newposition(2,ii));
                arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);newposition(2,ii)] [newposition(2,ii)+1;arcposition(2,newarcpos)+1] newarcmeshnum+arcposition(:,newarcpos+1:end)];
            end
            arcinterval=[arcinterval0(1:newarcpos) newarcintervaltime(2,ii) arcinterval(newarcpos+1:end)];
            y=[y0(:,1:newposition(1,ii)-1) ynew y(:,newposition(2,ii):end)];

            xnewr=x0(newposition(2,ii):arcposition0(2,newarcpos));
            xnewm=(xnewm-xnewm(1))/(xnewm(end)-xnewm(1))+newarcpos-1;
            xnewr=(xnewr-xnewr(1))/(xnewr(end)-xnewr(1))+newarcpos;
            if newarcpos<numarc
                x=[x0(:,1:newposition(1,ii)-1) xnewl xnewm xnewr 1+x(:,arcposition0(1,newarcpos+1):end)];
            else
                x=[x0(:,1:newposition(1,ii)-1) xnewl xnewm xnewr];
            end
            tgt=[];
        case 1
            newarcarg=[arcarg0(1:newarcpos) newarcid(ii) newarcarg(newarcpos+1:end)];
            newarcargindex=[arcargindex0(1:newarcpos) NaN newarcargindex(newarcpos+1:end)];
            xnewr=[];
            if newarcmeshnum>1
                xnewm=linspace(0,1,newarcmeshnum);
                if newposition(2,ii)==newposition(1,ii)
                    arcposition=[arcposition0(:,1:newarcpos) [arcposition0(2,newarcpos)+1;arcposition0(2,newarcpos)+newarcmeshnum] newarcmeshnum+arcposition(:,newarcpos+1:end)];
                else
                    arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);arcposition0(2,newarcpos)-1] [arcposition0(2,newarcpos);arcposition0(2,newarcpos)+newarcmeshnum-1] newarcmeshnum-1+arcposition(:,newarcpos+1:end)];
                end
            else
                xnewm=x0(newposition(1,ii):newposition(2,ii));
                arcposition=[arcposition0(:,1:newarcpos-1) [arcposition0(1,newarcpos);newposition(1,ii)] [newposition(1,ii)+1;arcposition(2,newarcpos)+1] newarcmeshnum+arcposition(:,newarcpos+1:end)];
            end
            arcinterval=[arcinterval0(1:newarcpos) newarcintervaltime(1,ii) arcinterval(newarcpos+1:end)];
            y=[y0(:,1:newposition(1,ii)) ynew y(:,newposition(2,ii)+1:end)];

            xnewl=x0(arcposition0(1,newarcpos):newposition(1,ii));
            xnewl=(xnewl-xnewl(1))/(xnewl(end)-xnewl(1))-1+newarcpos;
            xnewm=(xnewm-xnewm(1))/(xnewm(end)-xnewm(1))+newarcpos;
            if newarcpos>1
                x=[x0(:,1:arcposition0(2,newarcpos-1)) xnewl xnewm xnewr 1+x(:,newposition(2,ii)+1:end)];
            else
                x=[xnewl xnewm xnewr 1+x(:,newposition(2,ii)+1:end)];
            end
            tgt=[];
    end
end

sol.x=x;
sol.y=y;
sol.arcarg=newarcarg;
sol.arcposition=arcposition;
sol.arcinterval=arcinterval;
if ~isempty(tgt)
    sol.solverinfo.coeff=coeff;
    sol.solverinfo.tmesh=x;
    sol.solverinfo.tangent=tgt;
    sol.solverinfo.parametercoord=parametercoord;
    sol.solverinfo.meshvalcoord=meshvalcoord;
end