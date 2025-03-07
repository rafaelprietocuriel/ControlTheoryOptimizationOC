function ocTrj=octrajectoryprimitive(ocObj,X0,T,arcarg,opt)
%
% TRANSVERSALITYCONDITION

if isempty(ocObj)
    ocTrj=[];
    return
end

if nargin==3
    arcarg=0;
    opt=defaultocoptions;
end
if isempty(arcarg)
    arcarg=0;
end
if nargin==4
    opt=defaultocoptions;
end
if isempty(opt)
    opt=defaultocoptions;
end
TrivialArcMeshNum=getocoptions(opt,'GENERAL','TrivialArcMeshNum'); % number of initial time grid (repeated entries of equilibrium)

if T==0
    LT=transversalitycondition(ocObj,X0);
    ocTrj.y=repmat([X0;LT],1,4);
    ocTrj.arcarg=arcarg;
    ocTrj.x0=T;
    ocTrj.arcinterval=[0 0];
    ocTrj.timehorizon=0;
    ocTrj.modelname=modelname(ocObj);
    ocTrj.modelparameter=parametervalue(ocObj);
else
end

ocTrj.x=linspace(0,1,TrivialArcMeshNum);
ocTrj.y=repmat(ocTrj.y(:,1),1,TrivialArcMeshNum);
ocTrj.arcposition=[1;TrivialArcMeshNum];

ocTrj=octrajectory(ocTrj);
