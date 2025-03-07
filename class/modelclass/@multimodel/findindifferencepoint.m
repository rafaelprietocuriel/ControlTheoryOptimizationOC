function [indiffpt,o,idx1,idx2]=findindifferencepoint(mmObj,idx1,idx2,varargin)
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
if isempty(mmObj)
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
PlotResult=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
contRes=contresult(mmObj);

if nargin==2
    contSol1=contsolution(contRes{idx1});
    v1=zeros(1,length(contSol1));
    objval1=zeros(1,length(contSol1));
    if isfield(contSol1(1).solverinfo,'objectivevaluecoord')
        objectivevaluecoord=contSol1(1).solverinfo.objectivevaluecoord;
    else
        objectivevaluecoord=[];
    end
    arctimeindex=[];
    switch contSol1(1).solverinfo.continuationtype
        case 1
            arctimeindex=2;
            continuationindex=[];
    end
    if isfield(contSol1(1).solverinfo,'continuationindex')
        continuationindex=contSol1(1).solverinfo.continuationindex;
    else
        continuationindex=[];
    end
    for ii=1:length(contSol1)
        if ~isempty(continuationindex)
            v1(ii)=contSol1(ii).modelparameter{2}(continuationindex{1});
        elseif ~isempty(arctimeindex)
            v1(ii)=contSol1(ii).arcinterval(arctimeindex);
        end
        if ~isempty(objectivevaluecoord)
            objval1(ii)=contSol1(ii).y(objectivevaluecoord,end);
        end
    end
    startpt=0;
    v2=[];
    objval2=[];
elseif nargin==3
    contSol1=contsolution(contRes{idx1});
    contSol2=contsolution(contRes{idx2});
    v1=zeros(1,length(contSol1));
    objval1=zeros(1,length(contSol1));
    if isfield(contSol1(1).solverinfo,'objectivevaluecoord')
        objectivevaluecoord=contSol1(1).solverinfo.objectivevaluecoord;
    else
        objectivevaluecoord=[];
    end
    if isfield(contSol1(1).solverinfo,'continuationindex')
        continuationindex=contSol1(1).solverinfo.continuationindex;
        startpt=0;
    else
        continuationindex=[];
    end
    for ii=1:length(contSol1)
        if ~isempty(continuationindex)
            if iscell(continuationindex)
                v1(ii)=contSol1(ii).modelparameter{1}(continuationindex{1});
            else
                switch continuationindex
                    case 2
                        v1(ii)=contSol1(ii).solverinfo.initialtime4continuation+contSol1(ii).solverinfo.continuationparameter*contSol1(ii).solverinfo.continuationvector;
                end
            end
        end
        if ~isempty(objectivevaluecoord)
            objval1(ii)=contSol1(ii).y(objectivevaluecoord,end);
        end
    end

    %%%%
    v2=zeros(1,length(contSol2));
    objval2=zeros(1,length(contSol2));
    if isfield(contSol2(1).solverinfo,'objectivevaluecoord')
        objectivevaluecoord=contSol2(1).solverinfo.objectivevaluecoord;
    else
        objectivevaluecoord=[];
    end
    if isfield(contSol2(1).solverinfo,'continuationindex')
        continuationindex=contSol2(1).solverinfo.continuationindex;
    else
        continuationindex=[];
    end
    for ii=1:length(contSol2)
        if ~isempty(continuationindex)
            if iscell(continuationindex)
                v2(ii)=contSol2(ii).modelparameter{1}(continuationindex{1});
            else
                switch continuationindex
                    case 2
                        v2(ii)=contSol2(ii).solverinfo.initialtime4continuation+contSol2(ii).solverinfo.continuationparameter*contSol2(ii).solverinfo.continuationvector;
                end
            end
        end
        if ~isempty(objectivevaluecoord)
            objval2(ii)=contSol2(ii).y(objectivevaluecoord,end);
        end
    end
else
    return
end
if isempty(v2)
    [v0,o,idx1,idx2] = intersections(v1,objval1);
else
    [v0,o,idx1,idx2] = intersections(v1,objval1,v2,objval2);
end

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
    indiffpt(:,ii)=v0(ii);
    if PlotResult
        if isempty(v2)
            h=plot(v1,objval1,[v0(ii) v0(ii)],[min(objval1) max(objval1)]);
        else
            h=plot(v1,objval1,v2,objval2,[v0(ii) v0(ii)],[min(objval1) max(objval1)]);
        end
        if ii==1
            hold on
        elseif ii==length(v0)
            hold off
        end
        set(h,'Tag',['IPT' num2str(ii)])
        figure(gcf)
    end
end