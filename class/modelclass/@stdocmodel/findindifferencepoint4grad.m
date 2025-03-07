function [indiffval,o,idx1,idx2]=findindifferencepoint4grad(ocObj,contidx1,contidx2,varargin)
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

o=[];
opt=[];
if isempty(ocObj)
    return
end

if nargin>=4
    opt=varargin{1};
end

if isempty(opt)
    opt=defaultocoptions;
end
PlotResult=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');

if ~isnumeric(contidx1)
    return
end
contRes=contresult(ocObj);
contSol=contsolution(contRes{contidx1});
ocTrj=extremalsolution(ocObj);
contclass1=contclass(ocTrj{contidx1});
conttype1=conttype(ocTrj{contidx1});

if ~isocgradtrajectory(ocTrj{contidx1})
    return
end
c=contRes{contidx1}.ContinuationInformation;
solidx=length(contSol);
for ii=1:length(c)
    if strcmp(c(ii).label,'OVJ')
        solidx=c(ii).index;
    end
end
objval1=zeros(1,solidx);
contdata1=zeros(1,solidx);

for ii=1:solidx
    objval1(ii)=contSol(ii).extremal.objectivevalue;
    switch contclass1
        case 'extremalgrad4ft'
            if strcmp(conttype1,'time')
                contdata1(ii)=contSol(ii).extremal.timehorizon;
            elseif strcmp(conttype1,'initialstate')
                contdata1(ii)=contSol(ii).extremal.y(1,contSol(ii).extremal.initialstatecoordinate(1));
            elseif strcmp(conttype1,'parameter')
                contdata1(ii)=contSol(ii).extremal.modelparameter(contSol(ii).extremal.freeparameterindex);
            else
                return
            end
        case 'extremalgrad2ep'
            if strcmp(conttype1,'time')
                contdata1(ii)=contSol(ii).extremal.timehorizon;
            elseif strcmp(conttype1,'state')
                contdata1(ii)=contSol(ii).extremal.y(2,1);
            elseif strcmp(conttype1,'parameter')
                contdata1(ii)=contSol(ii).extremal.modelparameter(contSol(ii).extremal.freeparameterindex);
            else
                return
            end
            %ocgTrj=ocgradtrajectory(contSol(ii).extremal);
            %H=hamiltonian(ocgTrj);
            %objval1(ii)=H(1);

    end
end
if ~isnumeric(contidx2)
    return
end
contclass2=contclass(ocTrj{contidx2});
conttype2=conttype(ocTrj{contidx2});

if ~isocgradtrajectory(ocTrj{contidx1}) || ~strcmp(contclass1,contclass2) || ~strcmp(conttype1,conttype2)
    return
end
c=contRes{contidx2}.ContinuationInformation;
contSol=contsolution(contRes{contidx2});
solidx=length(contSol);
for ii=1:length(c)
    if strcmp(c(ii).label,'OVJ')
        solidx=c(ii).index;
    end
end
objval2=zeros(1,solidx);
contdata2=zeros(1,solidx);
for ii=1:solidx
    objval2(ii)=contSol(ii).extremal.objectivevalue;
    switch contclass2
        case 'extremalgrad4ft'
            if strcmp(conttype2,'time')
                contdata2(ii)=contSol(ii).extremal.timehorizon;
            elseif strcmp(conttype2,'initialstate')
                contdata2(ii)=contSol(ii).extremal.y(1,contSol(ii).extremal.initialstatecoordinate(1));
            elseif strcmp(conttype1,'parameter')
                contdata2(ii)=contSol(ii).extremal.modelparameter(contSol(ii).extremal.freeparameterindex);
            else
                return
            end
        case 'extremalgrad2ep'
            if strcmp(conttype2,'time')
                contdata2(ii)=contSol(ii).extremal.timehorizon;
            elseif strcmp(conttype2,'state')
                contdata2(ii)=contSol(ii).extremal.y(2,1);
            elseif strcmp(conttype2,'parameter')
                contdata2(ii)=contSol(ii).extremal.modelparameter(contSol(ii).extremal.freeparameterindex);
            else
                return
            end
            %ocgTrj=ocgradtrajectory(contSol(ii).extremal);
            %H=hamiltonian(ocgTrj);
            %objval2(ii)=H(1);
    end

end
[indiffval,o,idx1,idx2]=intersections(contdata1,objval1,contdata2,objval2);

 % remove spurious solutions 
 remidx=zeros(1,length(indiffval));
 for ii=1:length(indiffval)
     if abs(idx1(ii)-length(contdata1))<1
         remidx(ii)=ii;
     elseif abs(idx2(ii)-length(contdata2))<1
         remidx(ii)=ii;
     end
 end
remidx=find(remidx);
if ~isempty(remidx)
    indiffval(remidx)=[];
    o(remidx)=[];
    idx1(remidx)=[];
    idx2(remidx)=[];
end

if isempty(indiffval)
    if PlotResult
        plot(contdata1,objval1,contdata2,objval2)
        figure(gcf)
    end
    ocmatmsg('No indifference point detected.\n')
    return
end
for ii=1:length(indiffval)
    if PlotResult
        h=plot(contdata1,objval1,contdata2,objval2,[indiffval(ii) indiffval(ii)],[min([objval1 objval2]) max([objval1 objval2])]);
        if ii==1
            hold on
        elseif ii==length(indiffval)
            hold off
        end
        set(h,'Tag',['IPT' num2str(ii)])
        figure(gcf)
    end
end
