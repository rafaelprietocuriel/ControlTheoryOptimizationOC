function [indiffpt,o,idx1,idx2]=findindifferencepoint(ocObj,idx1,idx2,varargin)
%
% FINDINDIFFERENCEPOINT find indifference threshold.
%
% FINDINDIFFERENCEPOINT(OCOBJ,IDX1,IDX2) IDX1 and IDX2 are the indices of
% continuation results stored in OCOBJ. To detect an indifference threshold
% the Hamiltonian is evaluated at the initial points of the continuation
% process (slice manifold) and intersected. The slice manifold is only
% considered from the first point to the point where it bends back, if it
% bends back at all.
%
% FINDINDIFFERENCEPOINT(OCOBJ,SLMF1,SLMF2) SLMF1 and SLMF2 are slice
% manifolds.
%
% FINDINDIFFERENCEPOINT(OCOBJ,IDX1,IDX2,OPT)
%   OPT.GENERAL.AdmissibleTolerance defines the tolerance for the
%           collinearity test of the slice manifolds.
%   OPT.OCCONTARG.PlotCont='on'/'off' if set 'on' the result is shown
%           graphically
%
% INDIFFPT=FINDINDIFFERENCEPOINT(...) if the Hamiltonian functions
% intersect the (state) values of the intersection point INDIFFPT is
% returned otherwise it is empty.
%
% [INDIFFPT O]=FINDINDIFFERENCEPOINT(...) O is the corresponding objective
% value or empty.

indiffpt=[];
o=[];
opt=[];
contname='';
if isempty(ocObj)
    return
end

if nargin>=4
    opt=varargin{1};
end

if nargin>=5
    contname=varargin{2};
end

if isempty(opt)
    opt=defaultocoptions;
end
tol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
PlotResult=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');

contRes=contResult(ocObj);
if length(contRes)<max(idx1,idx2)
    return
end

for ii=1:order(ocObj)
    contRes{ii}=contresult(ocObj(ii));
end
if ismmultipath(ocTrj1{1})
    contRes=contResult(ocObj);
    contSol=contsolution(contRes{idx1});
    for ii=1:length(contSol)
        switch contSol(ii).solverinfo.continuationtype
            case 0
                v1(ii)=contSol(ii).solverinfo.continuationparameter;
            case 1
                v1(ii)=contSol(ii).solverinfo.continuationparameter;
            case 2
                v1(:,ii)=contSol(ii).solverinfo.initialparameter+contSol(ii).solverinfo.continuationparameter*contSol(ii).solverinfo.continuationvector;
                objval1(ii)=contSol(ii).y(contSol(ii).solverinfo.objectivevaluecoord,end);
                refvec=contSol(ii).solverinfo.continuationparameter;
                startpt=contSol(ii).solverinfo.initialparameter;
        end

    end
    v1=v1(1,:);
else
    state1=state(ocObj,idx1);
    if size(state1,1)>1
        nrm1=sqrt(sum(state1.^2));
    else
        nrm1=abs(state1);
    end
    [dum idx]=min(nrm1);
    startpt1=state1(:,idx);
    [dum idx]=max(nrm1);
    endpt1=state1(:,idx);
end
[v0,o,idx1,idx2] = intersections(v1,objval1);


if isempty(v0)
    if PlotResult
        if ~isempty(sliceman2)
            plot(v1,objval1,v2,objval2)
        else
            plot(v1,objval1)
        end
        figure(gcf)
    end
    ocmatmsg('No indifference point detected.\n')
    return
end

indiffpt=zeros(size(startpt,1),length(v0));
for ii=1:length(v0)
    indiffpt(:,ii)=v0(ii)*refvec;
    if PlotResult
        h=plot(v1,objval1,[v0(ii) v0(ii)],[min(objval1) max(objval1)]);
        if ii==1
            hold on
        elseif ii==length(v0)
            hold off
        end
        set(h,'Tag',['IPT' num2str(ii)])
        figure(gcf)
    end
end