function [indiffpt,o,idx1,idx2]=findindifferencepoint(docObj,contidx1,contidx2,varargin)
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
% [INDIFFPT O]=FINDINDIFFERENCEPOINT(...) O is the corresponding obective
% value or empty.

indiffpt=[];
o=[];
opt=[];
lightgray=repmat(0.75,1,3);
black=[0 0 0];
red=[1 0 0];
if isempty(docObj)
    return
end
if nargin==2
    contidx2=[];
end
if nargin>=4
    opt=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end
tol=getocoptions(opt,'GENERAL','AdmissibleTolerance');
PlotResult=strcmpi(getocoptions(opt,'OCCONTARG','PlotCont'),'on');
statecrd=statecoord(docObj);
statenum=length(statecrd);
sliceman1=slicemanifold(docObj,contidx1);
if ~isempty(sliceman1)
    ContinuationSolution1=docObj.Result.Continuation{contidx1}.ContinuationSolution;
    docFP1=docObj.Result.Continuation{contidx1}.LimitSet;
    if isempty(contidx2)
        sliceman2=occurve([]);
    else
        ContinuationSolution2=docObj.Result.Continuation{contidx2}.ContinuationSolution;
        docFP2=docObj.Result.Continuation{contidx2}.LimitSet;
        sliceman2=slicemanifold(docObj,contidx2);
    end
    if ~isempty(sliceman2)
        state1=state(docObj,sliceman1);
        state2=state(docObj,sliceman2);
        objval1=sliceman1.y(end,:);
        objval2=sliceman2.y(end,:);

        % test if slice manifolds are collinear in the state space
        if ~(rank([diff(state1(:,[1 end]),[],2) diff(state2(:,[1 end]),[],2)],tol)==1)
            ocmaterror('Directions of slice manifolds are not colinear.')
        end
        statept=[state1 state2];
        if size(statept,1)>1
            normstate=sum(statept.^2);
            normstate1=sqrt(sum(state1.^2));
            normstate2=sqrt(sum(state2.^2));
        else
            normstate=statept.^2;
            normstate1=sqrt(state1.^2);
            normstate2=sqrt(state2.^2);
        end
        [dum minidx]=min(normstate);
        [dum maxidx]=max(normstate);
        startpt=statept(:,minidx);
        endpt=statept(:,maxidx);
        refvec=endpt-startpt;
        %         if PlotResult
        %             clf
        %             h=plot(normstate1,objval1,normstate2,objval2);
        %             set(h(1),'Color',lightgray);
        %             set(h(2),'Color',black);
        %             hold on
        %             minobjval=min([objval1,objval2]);
        %             maxobjval=max([objval1,objval2]);
        %         end

        v1=pinv(refvec)*(state1-startpt(:,ones(1,size(state1,2))));
        v2=pinv(refvec)*(state2-startpt(:,ones(1,size(state2,2))));
        [v0,o,idx1,idx2]=intersections(v1,objval1,v2,objval2);

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
    end
else
    if ~exist('refvec','var')
        startpt=startpt1;
        refvec=endpt1-startpt;
    end
    v1=continuationparameter(sliceman1);
    [v0,o,idx1,idx2] = intersections(v1,objval1);
end
indiffpt=zeros(size(startpt,1),length(v0));
for ii=1:length(v0)
    indiffpt(:,ii)=startpt+v0(ii)*refvec;
    if PlotResult
        if ~isempty(sliceman2)
            h=plot(v1,objval1,v2,objval2,[v0(ii) v0(ii)],[min([objval1 objval2]) max([objval1 objval2])]);
        else
            h=plot(v1,objval1,[v0(ii) v0(ii)],[min(objval1) max(objval1)]);
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
