function [indiffpt o]=findindifferencedistribution(ppdeObj,sliceman1,sliceman2,varargin)
%
% FINDINDIFFERENCEDISTRIBUTION find indifference threshold.
%
% FINDINDIFFERENCEDISTRIBUTION(OCOBJ,IDX1,IDX2) IDX1 and IDX2 are the indices of
% continuation results stored in OCOBJ. To detect an indifference threshold
% the Hamiltonian is evaluated at the initial points of the continuation
% process (slice manifold) and intersected. The slice manifold is only
% considered from the first point to the point where it bends back, if it
% bends back at all.
%
% FINDINDIFFERENCEDISTRIBUTION(OCOBJ,IDX1,IDX2,OPT) 
%   OPT.GENERAL.AdmissibleTolerance defines the tolerance for the
%           collinearity test of the slice manifolds.  
%   OPT.OCCONTARG.PlotCont='on'/'off' if set 'on' the result is shown
%           graphically
%
% FINDINDIFFERENCEDISTRIBUTION(OCOBJ,IDX1,IDX2,OPT,SPEC1,SPEC2) the argument
%       SPEC1= ''/'r': if set to 'r' the slice manifold is considered in
%               reverse order. This argument is useful if the threshold
%               point separates two solutions converging to the same limit
%               set. In that case it suffices to compute one slice manifold
%               (with at least two limit points) IDX and call
%               FINDINDIFFERENCEDISTRIBUTION(OCOBJ,IDX,IDX,[],'r')
%
% INDIFFPT=FINDINDIFFERENCEDISTRIBUTION(...) if the Hamiltonian functions
% intersect the (state) values of the intersection point INDIFFPT is
% returned otherwise it is empty.
%
% [INDIFFPT O]=FINDINDIFFERENCEDISTRIBUTION(...) O is the corresponding obective
% value or empty.

indiffpt=[];
opt=[];
spec1='';
spec2='';

if isempty(ppdeObj)
    return
end

if nargin>=4
    opt=varargin{1};
end
if nargin>=5
    spec1=varargin{2};
end
if nargin>=6
    spec2=varargin{3};
end

if isempty(opt)
    opt=defaultocoptions;
end
tol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
PlotResult=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
contRes=contresult(ppdeObj);
if isnumeric(sliceman1)
    counter=1;
    len=length(contRes{sliceman1}.ContinuationSolution);
    objval1=zeros(1,len);
    contpar1=objval1;
    for ii=1:len
        try
            ppdeTrj=ppdeasymptotic(ppdetrajectory(contRes{sliceman1}.ContinuationSolution(ii)),contRes{sliceman1}.LimitSet);
            objval1(counter)=objectivevalue(ppdeObj,ppdeTrj);
            contpar1(counter)=ppdeTrj.freeparameter(end);
            state(ppdeObj,ppdeTrj);
            state1(:,counter)=ans(:,1);
            counter=counter+1;
        end
    end
elseif iscell(sliceman1) && isoccurve(sliceman1{1})
    return
end

if isnumeric(sliceman2)
    counter=1;
    objval2=zeros(1,length(contRes{sliceman2}.ContinuationSolution));
    contpar2=objval2;
    for ii=1:length(contRes{sliceman2}.ContinuationSolution)
        try
            ppdeTrj=ppdeasymptotic(ppdetrajectory(contRes{sliceman2}.ContinuationSolution(ii)),contRes{sliceman2}.LimitSet);
            objval2(counter)=objectivevalue(ppdeObj,ppdeTrj);
            contpar2(counter)=ppdeTrj.freeparameter(end);
            state(ppdeObj,ppdeTrj);
            state2(:,counter)=ans(:,1);
            counter=counter+1;
        end
    end
elseif iscell(sliceman2) && isoccurve(sliceman2{1})
    return
end


doubleentryidx1=find(diff(contpar1)==0);
contpar1(doubleentryidx1)=[];
contpar1idx=1:numel(contpar1);
doubleentryidx2=find(diff(contpar2)==0);
contpar2(doubleentryidx2)=[];
contpar2idx=1:numel(contpar2);

% points have to be considered in the reversed order
if strcmp(spec1,'r')
    contpar1=fliplr(contpar1);
    contpar1idx=fliplr(contpar1idx);
end

% points have to be considered in the reversed order
if strcmp(spec2,'r')
    contpar2=fliplr(contpar2);
    contpar2idx=fliplr(contpar2idx);
end

if numel(contpar1)<=3 || numel(contpar2)<=3
    ocmatmsg('Slice manifold consits of less than four points.\n')
    return
end


% find index where direction changes its sign
idx1=find(diff(sign(diff(contpar1)))&abs(contpar1(1:end-2))>1e-4,1);
idx2=find(diff(sign(diff(contpar2)))&abs(contpar2(1:end-2))>1e-4,1);
if isempty(idx1)
    idx1=numel(contpar1idx);
end
objval1=objval1(contpar1idx(1:idx1));
state1=state1(:,contpar1idx(1:idx1));

if isempty(idx2)
    idx2=numel(contpar2idx);
end
objval2=objval2(contpar2idx(1:idx2));
state2=state2(:,contpar2idx(1:idx2));
% test if slice manifolds are collinear in the state space
if ~(rank([diff(state1(:,[1 end]),[],2) diff(state2(:,[1 end]),[],2)],tol)==1)
    ocmaterror('Directions of slice manifolds are not colinear.')
end
startpt=state1(:,1);
endpt=state2(:,1);
refvec=endpt-startpt;

v1=pinv(refvec)*(state1-startpt(:,ones(1,idx1)));
v2=pinv(refvec)*(state2-startpt(:,ones(1,idx2)));

v=union(v1,v2);

objval1int=interp1(v1,objval1,v,'linear');
objval2int=interp1(v2,objval2,v,'linear');

diffobjval=objval2int-objval1int;
removeidx=isnan(diffobjval);
diffobjval(removeidx)=[];
v(removeidx)=[];
% find changes in the sign
idx=find(diff(sign(diffobjval)),1);

if isempty(idx)
    ocmatmsg('No indifference point detected.')
    return
end

o=diffobjval(idx:idx+1);
do=diff(o);
v=v(idx:idx+1);
if abs(do)<1e-10
    v0=(v(1)+v(2))/2;
else
    v0=v(1)-o(1)/(o(2)-o(1))*(v(2)-v(1));
end
indiffpt=startpt+v0*refvec;
if PlotResult
    plot(v1,objval1,v2,objval2,[v0 v0],[min([objval1int objval2int]) max([objval1int objval2int])])
    figure(gcf)
end
