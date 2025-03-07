function [hocTrj,info]=hybridprimitive(ocObj,X0,T,jumparg,arcarg,opt)
%
% HYBRIDPRIMITIVE returns an initial constant hybridoctrajectory
%
% HOCTRJ=HYBRIDPRIMITIVE(OCOBJ,X0,T,JUMPARG) information is taken from the
% impulse model OCOBJ. An hybridoctrajectory HOCTRJ starting at X0 being
% constant on the time interval T is returned. The transversality condition
% is calculated for S(X0,T(end)). In general this is only correct for T=0.
% The jump argument JUMPARG provides inforamtion about jumps at the
% beginning or end: 
%   [0 0]       no jumps
%   [1 0]       jump at the beginning
%   [0 1]       jump at the end
%   [1 1]       jump at the beginning and the end
%
% HYBRIDPRIMITIVE(OCOBJ,X0,T,JUMPARG,ARCARG) if the impulse model OCOBJ
% includes a usual control ARCARG provides information about active and
% inactive constraints.
%
% HYBRIDPRIMITIVE(OCOBJ,X0,T,JUMPARG,ARCARG,OPT) the option
% 'TrivialArcMeshNum' in the 'GENERAL' field returns the number for the time grid.  
%
% [HOCTRJ,INFO]=HYBRIDPRIMITIVE(OCOBJ,X0,T,JUMPARG,...) INFO returns
% information about the 

if isempty(ocObj)
    hocTrj=[];
    return
end

if nargin==4
    arcarg=0;
    opt=defaultocoptions;
end
if isempty(arcarg)
    arcarg=0;
end
if nargin==5
    opt=defaultocoptions;
end
if isempty(opt)
    opt=defaultocoptions;
end
if length(T)==1
    T=[0 T];
end
info=[];
TrivialArcMeshNum=getocoptions(opt,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)

X0=X0(:);
LT=repmat(NaN,size(X0));%transversalitycondition(ocObj,X0);
hocTrj.y=repmat([X0;LT],1,4);
hocTrj.jumparg=jumparg;
hocTrj.arcarg=arcarg;
hocTrj.x0=0;
hocTrj.arcinterval=T;
hocTrj.timehorizon=T(end);
hocTrj.modelname=modelname(ocObj);
hocTrj.modelparameter=parametervalue(ocObj);

LT=transversalitycondition(ocObj,hocTrj);
hocTrj.y=repmat([X0;LT],1,4);

if any(jumparg)
    modeln=modelname(ocObj);
    ics.jump=str2func([modeln 'EventBC']);
    ics.transversality=str2func([modeln 'TransversalityBC']);
    [depvar fval fflag info]=fsolve(@zerofunc,hocTrj.y(:),opt.EQ,[0 T],hocTrj.modelparameter,X0,hocTrj.arcarg,hocTrj.jumparg,ics);
    if ~fflag==1
        ocmatmsg('Convergence to a correct solution is not guaranteed.\n')
    end
    hocTrj.y=reshape(depvar,size(hocTrj.y));
end
hocTrj.x=[0 linspace(0,1,TrivialArcMeshNum) 1];
hocTrj.y=[hocTrj.y(:,1) repmat(hocTrj.y(:,2),1,TrivialArcMeshNum) hocTrj.y(:,end)];
hocTrj.arcposition=[2;TrivialArcMeshNum+1];

hocTrj=hybridoctrajectory(hocTrj);

    function res=zerofunc(depvar,timepoints,par,X0,arcarg,jumparg,ics)

        depvar=reshape(depvar,[],4);
        res=depvar(1:length(X0),1)-X0;
        for ii=1:2
            depvarLR=depvar(:,(2*ii-1):2*ii); % variable representing the left and right limits
            if ii==1
                % initial impulse or continuity
                res=[res; ...
                    ics.jump(timepoints(ii),depvarLR,par,arcarg,jumparg(ii))];
                depvarI=depvar(:,2*ii:(2*ii+1)); % variable representing the points on the interval between two impulses
                res=[res;diff(depvarI,[],2)];
            elseif ii==2
                % end impulse or continuity
                res=[res; ...
                    ics.jump(timepoints(ii),depvarLR,par,arcarg,jumparg(ii))];
                % transversality condition
                res=[res;ics.transversality(timepoints(ii),depvarLR,par,arcarg,jumparg(ii))];
            end
        end
    end
end