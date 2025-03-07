function [indiffpt,o,solidx1,solidx2]=findstaticindifferencepoint(ocObj,idx1,idx2,varargin)
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

continuationtype='';
coordinate=[];
parameter='';
option='';


if isempty(ocObj)
    return
end
optionidx=find(strcmpi(varargin,'option'));
coordinateidx=find(strcmpi(varargin,'coordinate'));
continuationtypeidx=find(strcmpi(varargin,'continuationtype'));

if ~isempty(optionidx)
    option=varargin{optionidx+1};
end
if ~isempty(coordinateidx)
    coordinate=varargin{coordinateidx+1};
end
if ~isempty(continuationtypeidx)
    continuationtype=varargin{continuationtypeidx+1};
end
if strcmp(continuationtype,'parameter')
    paridx=parameterindex(ocObj,varargin{continuationtypeidx+2});
end
if isempty(option)
    opt=defaultocoptions;
end
PlotResult=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');


contRes1=contresult(ocObj,idx1);
contSol1=contsolution(contRes1{1});
contRes2=contresult(ocObj,idx2);
contSol2=contsolution(contRes2{1});

optval1=zeros(1,length(contSol1));
contval1=zeros(1,length(contSol1));
optval2=zeros(1,length(contSol2));
contval2=zeros(1,length(contSol2));

for ii=1:length(contSol1)
    if ~isempty(findstr(contSol1(ii).solverinfo.conttype,'indifferencesolution'))
    elseif ~isempty(findstr(contSol1(ii).solverinfo.conttype,'4ft'))
        ocTrj1=octrajectory(contSol1(ii));
        depvar1=dependentvar(ocTrj1);
        optval1(ii)=depvar1(coordinate,end);
    else
    end
    switch continuationtype
        case 'parameter'
            par=modelparameter(ocTrj1);
            contval1(ii)=par(paridx);
    end
end

for ii=1:length(contSol2)
    if ~isempty(findstr(contSol2(ii).solverinfo.conttype,'indifferencesolution'))
    elseif ~isempty(findstr(contSol2(ii).solverinfo.conttype,'4ft'))
        ocTrj2=octrajectory(contSol2(ii));
        depvar2=dependentvar(ocTrj2);
        optval2(ii)=depvar2(coordinate,end);
    else
    end
    switch continuationtype
        case 'parameter'
            par=modelparameter(ocTrj2);
            contval2(ii)=par(paridx);
    end
end

[indiffpt,o,solidx1,solidx2]=intersections(contval1,optval1,contval2,optval2);
remidx=[find(isnan(solidx1));find(isnan(solidx2))];
indiffpt(remidx)=[];
solidx1(remidx)=[];
solidx2(remidx)=[];

if isempty(indiffpt)
    if PlotResult
        if ~isempty(idx2)
            plot(contval1,optval1,contval2,optval2)
        else
            plot(contval1,optval1)
        end
        figure(gcf)
    end
    ocmatmsg('No indifference point detected.\n')
    return
end

for ii=1:length(indiffpt)
    if PlotResult
        if ~isempty(idx2)
            h=plot(contval1,optval1,contval2,optval2,[indiffpt(ii) indiffpt(ii)],[min([optval1 optval2]) max([optval1 optval2])]);
        else
            h=plot(contval1,optval1,[indiffpt(ii) indiffpt(ii)],[min(optval1) max(optval1)]);
        end
        if ii==1
            hold on
        elseif ii==length(indiffpt)
            hold off
        end
        set(h,'Tag',['IPT' num2str(ii)])
        figure(gcf)
    end
end