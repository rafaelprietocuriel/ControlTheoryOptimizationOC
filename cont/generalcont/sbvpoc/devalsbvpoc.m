function ynew=devalsbvpoc(coeff,tmesh,tnew,drearrfunc,domainddata)

global OCMATCONT

if nargin==3
    domainddata=OCMATCONT.DOMAINDDATA;
    drearrfunc=OCMATCONT.drearr;
    OCMATCONT.meshadaptationflag=1;
end

if nargin==4
    domainddata=OCMATCONT.DOMAINDDATA;
    OCMATCONT.meshadaptationflag=1;
end
[tnew sortidx]=sort(tnew);

% this function returns the values of the coefficients in terms of the y
% and z coordinates as explained in the bvpsuite manual.
[tmesh,y,z]=drearrfunc(tmesh,coeff);
diffmesh=diff(tmesh);

difftnew=diff(tnew);
arcindex=find(difftnew==0);
leftarcindexnew=[1 arcindex+1];
rightarcindexnew=[arcindex numel(tnew)];
numarcnew=numel(arcindex)+1;
numarc=OCMATCONT.HE.numarc;
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;

if sum(ismember(tnew,[tmesh(leftarcindex(2:numarc))]))-2*(numarcnew-1)
    error('At switching times the interpolated mesh points have to appear twice.')
end
if numarcnew>1 && ~all(ismember(tnew(leftarcindexnew(2:end)),tmesh(leftarcindex))) && ~all(ismember(tnew(rightarcindexnew(1:end-1)),tmesh(rightarcindex)))
    error('Switching times are inconsistent')
end
endpointflag=0;
for arcnew=1:numarcnew
    for ii=leftarcindexnew(arcnew):rightarcindexnew(arcnew)
        % handle the different possible cases for left and right limits.
        % The user has to assure that time points for the interpolated time
        % mesh at the switching times have to be doubled.
        % 'endpointflag' marks time points at the right end of the arc time
        % interval
        arc=find(tnew(ii)==tmesh(leftarcindex));
        if isempty(arc)
            if tnew(ii)==tmesh(rightarcindex(numarc))
                arc=numarc;
                endpointflag=1;
            else
                endpointflag=0;
                arc=find(tnew(ii)>tmesh(leftarcindex));
                arc=arc(end);
            end
        elseif arc>1
            if endpointflag
                endpointflag=0;
            else
                arc=arc(1)-1;
                endpointflag=1;
            end
        else
            endpointflag=0;
        end
        arcindex=OCMATCONT.HE.arcindex(arc);
        numcolscoord=domainddata(arcindex).numcolscoord;
        numae=domainddata(arcindex).numae;
        odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
        diffmesh_arc=diffmesh(leftarcindex(arc):rightarcindex(arc)-1);
        tmesh_arc=tmesh(leftarcindex(arc):rightarcindex(arc));
        y_arc=y(:,:,leftarcindex(arc)-arc+1:rightarcindex(arc)-arc);
        z_arc=z(:,numcolscoord,leftarcindex(arc)-arc+1:rightarcindex(arc)-arc);

        idx=find(tnew(ii)>=tmesh_arc);
        idx=idx(end);
        if ~endpointflag
            ynew(odecoord,ii)=evalcollocationpoly(idx,odecoord,tnew(ii),tmesh_arc,y_arc,z_arc,diffmesh_arc,numcolscoord,domainddata(arcindex).psi,1);
            if numae
                aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations
                ynew(aecoord,ii)=evalcollocationpoly(idx,aecoord,tnew(ii),tmesh_arc,y_arc,z_arc,diffmesh_arc,numcolscoord,domainddata(arcindex).psi0,0);
            end
        else
            % case where time points are lying on the right side of the arc
            % time interval. For mesh pointsin the interior of the arc time
            % interval the solutions have to be continuous and therefore
            % the left and right side limit coincide. At switching times the
            % left and right side limits may be different.
            ynew(odecoord,ii)=evalcollocationpoly(idx-1,odecoord,tnew(ii),tmesh_arc,y_arc,z_arc,diffmesh_arc,numcolscoord,domainddata(arcindex).psi,1);
            if numae
                aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations
                ynew(aecoord,ii)=evalcollocationpoly(idx-1,aecoord,tnew(ii),tmesh_arc,y_arc,z_arc,diffmesh_arc,numcolscoord,domainddata(arcindex).psi0,0);
            end
        end
    end
end
ynew=ynew(:,sortidx);