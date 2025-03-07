function [indifftime,o,idx1,idx2,hocTrj]=findindifferencetime(ocObj,sliceman1,sliceman2,varargin)
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

indifftime=[];
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

contRes=result(ocObj);
contSol=contRes.Continuation{sliceman1};
objval1=zeros(1,length(contSol.ContinuationSolution));
v1=objval1;
for ii=1:length(contSol.ContinuationSolution)
    hocTrj=hybridoctrajectory(contSol.ContinuationSolution(ii));
    tau=arcinterval(hocTrj);
    objval1(:,ii)=objectivevalue(ocObj,hocTrj);
    v1(ii)=tau(end);
end
contSol=contRes.Continuation{sliceman2};
objval2=zeros(1,length(contSol.ContinuationSolution));
v2=objval2;
for ii=1:length(contSol.ContinuationSolution)
    hocTrj=hybridoctrajectory(contSol.ContinuationSolution(ii));
    tau=arcinterval(hocTrj);
    objval2(:,ii)=objectivevalue(ocObj,hocTrj);
    v2(ii)=tau(end);
end
[indifftime,o,idx1,idx2] = intersections(v1,objval1,v2,objval2);


if isempty(indifftime)
    if PlotResult
        plot(v1,objval1,v2,objval2)
        figure(gcf)
    end
    ocmatmsg('No indifference point detected.\n')
    return
end

idx=[sliceman1 sliceman2];
indiffidx=round([idx1 idx2]);
hocTrj=cell(1,2);
for ii=1:length(indifftime)
    for jj=1:length(idx)
        contSol=contRes.Continuation{idx(jj)};
        if length(indifftime)==1
            hocTrj{jj}=hybridoctrajectory(contSol.ContinuationSolution(indiffidx(ii,jj)));
        end
    end
    if PlotResult
        plot(v1,objval1,v2,objval2,[indifftime(ii) indifftime(ii)],[min([objval1 objval2]) max([objval1 objval2])])
        if ii==1
            hold on
        elseif ii==length(indifftime)
            hold off
        end
        figure(gcf)
    end
end