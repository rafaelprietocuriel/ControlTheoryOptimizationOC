function [indiffpt,o,idx1,idx2]=findindifferencepoint(ocObj,sliceman1,sliceman2,varargin)
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
isfinitesol1=[];
isfinitesol2=[];
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

if isnumeric(sliceman1)
    ocTrj1=extremalsolution(ocObj,sliceman1,contname);
    if isoctrajectory(ocTrj1{1}) && isocasymptotic(ocTrj1{1})
        if statenum(ocObj)==1
            sliceman1=occurve(octrajectory(ocTrj1{1}));
        else
            sliceman1=occurve(slicemanifold(ocObj,sliceman1,contname));
        end
        isfinitesol1=0;
    elseif isa(ocTrj1{1},'ocgtrajectory') && isa(ocTrj1{1},'ocgasymptotic')
        contRes=contresult(ocObj);
        contSol1=contsolution(contRes{sliceman1});
        endman1=occurve(endpointmanifold(ocObj,sliceman1,contname));
        sliceman1=occurve(slicemanifold(ocObj,sliceman1,contname));
        isfinitesol1=0;
    elseif isoctrajectory(ocTrj1{1})
        contRes=contresult(ocObj);
        contSol1=contsolution(contRes{sliceman1});
        endman1=occurve(endpointmanifold(ocObj,sliceman1,contname));
        sliceman1=occurve(slicemanifold(ocObj,sliceman1,contname));
        isfinitesol1=1;
    elseif isocmultipath(ocTrj1{1})
        contRes=contresult(ocObj);
        contSol1=contsolution(contRes{sliceman1});
        mf=endpointmanifold(ocObj,sliceman1,contname);
        endman1=occurve(mf{1});
        mf=slicemanifold(ocObj,sliceman1,contname);
        sliceman1=occurve(mf{1});
        isfinitesol1=1;
    end
elseif iscell(sliceman1) && isoccurve(sliceman1{1})
    sliceman1=sliceman1{1};
elseif iscell(sliceman1) && isoctrajectory(sliceman1{1}) && isocasymptotic(sliceman1{1})
    sliceman1=occurve(octrajectory(sliceman1{1}));
elseif isocasymptotic(sliceman1)
    sliceman1=occurve(octrajectory(sliceman1));
end
if ~isoccurve(sliceman1)
    ocmaterror('Second argument does not yield a slice manifold.')
end

if ~isempty(sliceman2) && isnumeric(sliceman2)
    ocTrj2=extremalsolution(ocObj,sliceman2,contname);
    if isoctrajectory(ocTrj2{1}) && isocasymptotic(ocTrj2{1})
        if statenum(ocObj)==1
            sliceman2=occurve(octrajectory(ocTrj2{1}));
        else
            sliceman2=occurve(slicemanifold(ocObj,sliceman2,contname));
        end
        isfinitesol2=0;
    elseif isa(ocTrj2{1},'ocgtrajectory') && isa(ocTrj2{1},'ocgasymptotic')
        contRes=contresult(ocObj);
        contSol2=contsolution(contRes{sliceman2});
        endman2=occurve(endpointmanifold(ocObj,sliceman2,contname));
        sliceman2=occurve(slicemanifold(ocObj,sliceman2,contname));
        isfinitesol1=0;
    elseif isoctrajectory(ocTrj2{1})
        contRes=contresult(ocObj);
        contSol2=contsolution(contRes{sliceman2});
        endman2=occurve(endpointmanifold(ocObj,sliceman2,contname));
        sliceman2=occurve(slicemanifold(ocObj,sliceman2,contname));
        isfinitesol2=1;
    elseif isocmultipath(ocTrj2{1})
        contRes=contresult(ocObj);
        contSol2=contsolution(contRes{sliceman2});
        mf=endpointmanifold(ocObj,sliceman2,contname);
        endman2=occurve(mf{1});
        mf=slicemanifold(ocObj,sliceman2,contname);
        sliceman2=occurve(mf{1});
        isfinitesol2=1;
    end
elseif iscell(sliceman2) && isoccurve(sliceman2{1})
    sliceman2=sliceman2{1};
elseif iscell(sliceman2) && isoctrajectory(sliceman2{1})
    sliceman2=occurve(octrajectory(sliceman2{1}));
elseif isocasymptotic(sliceman2)
    sliceman2=occurve(octrajectory(sliceman2));
end

if ~isempty(sliceman2) && ~isoccurve(sliceman2)
    ocmatmsg('Third argument does not yield a slice manifold.')
    return
end

switch conttype(ocTrj1{1})
    case 'extremalp4ft'
        v1=zeros(1,length(contSol1));
        for ii=1:length(contSol1)
            v1(ii)=contSol1(ii).solverinfo.initialparametervaluevalue(1)+contSol1(ii).solverinfo.continuationparameter*contSol1(ii).solverinfo.continuationvector(1);
        end
        %v1=continuationparameter(sliceman1);
        startpt=v1(1);
        refvec=1;
    otherwise
        state1=state(ocObj,sliceman1,1);
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
if ~isfield(sliceman1.userinfo,'objectivevalue')
    if ~isfinitesol1
        r=double(discountrate(ocObj));
        objval1=hamiltonian(ocObj,sliceman1,1)/r;
    else
        r=double(discountrate(ocObj));
        if strcmp(conttype(ocTrj1{1}),'extremalt4ft')
            T=sliceman1.userinfo.continuationparameter;
            objval1=(hamiltonian(ocObj,sliceman1)-exp(-r*T).*hamiltonian(ocObj,endman1))/r+exp(-r*T).*salvagevalue(ocObj,endman1);
        elseif isocmultipath(ocTrj1{1})
            objval1=sliceman1.y(end,:);
        else
            T=timehorizon(ocTrj1{1});
            objval1=(hamiltonian(ocObj,sliceman1)-exp(-r*T)*hamiltonian(ocObj,endman1))/r+exp(-r*T)*salvagevalue(ocObj,endman1);
        end
    end
else
    objval1=sliceman1.userinfo.objectivevalue;
end

if ~isempty(sliceman2)
    if ~isfield(sliceman2.userinfo,'objectivevalue')
        if ~isfinitesol2
            r=double(discountrate(ocObj));
            objval2=hamiltonian(ocObj,sliceman2,1)/r;
        else
            r=double(discountrate(ocObj));
            if strcmp(conttype(ocTrj2{1}),'extremalt4ft')
                T=sliceman1.userinfo.continuationparameter;
                objval2=(hamiltonian(ocObj,sliceman2)-exp(-r*T).*hamiltonian(ocObj,endman2))/r+exp(-r*T)*salvagevalue(ocObj,endman2);
        elseif isocmultipath(ocTrj2{1})
            objval2=sliceman2.y(end,:);
            else
                T=timehorizon(ocTrj2{1});
                objval2=(hamiltonian(ocObj,sliceman2)-exp(-r*T)*hamiltonian(ocObj,endman2))/r+exp(-r*T)*salvagevalue(ocObj,endman2);
            end
        end
    else
        objval2=sliceman2.userinfo.objectivevalue;
    end
    switch conttype(ocTrj2{1})
        case 'extremalp4ft'
        v2=zeros(1,length(contSol2));
        for ii=1:length(contSol2)
            v2(ii)=contSol2(ii).solverinfo.initialparametervaluevalue(1)+contSol2(ii).solverinfo.continuationparameter*contSol2(ii).solverinfo.continuationvector(1);
        end
            %v2=continuationparameter(sliceman2);
        otherwise
            state2=state(ocObj,sliceman2,1);
            if size(state2,1)>1
                nrm2=sqrt(sum(state2.^2));
            else
                nrm2=abs(state2);
            end

            [dum idx]=min(nrm2);
            startpt2=state2(:,idx);
            [dum idx]=max(nrm2);
            endpt2=state2(:,idx);
            % test if slice manifolds are collinear in the state space
            if 0%~(rank([endpt1-startpt1 endpt2-startpt2],tol)==1)
                ocmatmsg('Directions of slice manifolds are not colinear.\n')
                return
            end
            refvec=endpt1-startpt1;
            startpt=startpt1;
            v1=pinv(refvec)*(state1-startpt(:,ones(1,size(state1,2))));
            v2=pinv(refvec)*(state2-startpt(:,ones(1,size(state2,2))));
    end
    [v0,o,idx1,idx2]=intersections(v1,objval1,v2,objval2);
else
    if ~exist('refvec','var')
        startpt=startpt1;
        refvec=endpt1-startpt;
    end
    v1=continuationparameter(sliceman1);
    [v0,o,idx1,idx2] = intersections(v1,objval1);
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

indiffpt=v0;
for ii=1:length(v0)
    %indiffpt(:,ii)=startpt+v0(ii)*refvec;
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