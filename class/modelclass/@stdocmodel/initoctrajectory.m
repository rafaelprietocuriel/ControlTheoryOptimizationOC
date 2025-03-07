function ocTrj=initoctrajectory(ocObj,x0,arcarg,varargin)

fixinitialstatecoord=[];
fixendstatecoord=[];
opt=[];

optionidx=find(strcmpi(varargin,'option'));
fixinitialstatecoordidx=find(strcmpi(varargin,'fixinitialstatecoord'));
fixendstatecoordidx=find(strcmpi(varargin,'fixendstatecoord'));

xcoord=statecoord(ocObj);
x0=x0(:);

if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if ~isempty(fixinitialstatecoordidx)
    fixinitialstatecoord=varargin{fixinitialstatecoordidx+1};
end
if ~isempty(fixendstatecoordidx)
    fixendstatecoord=varargin{fixendstatecoordidx+1};
end
if isempty(arcarg)
    arcarg=0;
end
if isempty(opt)
    opt=defaultocoptions;
end

n=getocoptions(opt,'GENERAL','TrivialArcMeshNum');

fixinitialcostatecoord=setdiff(xcoord,fixinitialstatecoord);
fixendcostatecoord=setdiff(xcoord,fixendstatecoord);


initocPt.y=x0;
initocPt.x=0;
initocPt.arcarg=arcarg;
initocPt.arcinterval=[0 0];
initocPt.arcposition=[1;1];

initialcostate(xcoord)=0;
if ~isempty(fixendcostatecoord)
    tmp=transversalitycondition(ocObj,octrajectory(initocPt));
    initialcostate(fixendcostatecoord)=tmp(fixendcostatecoord);
end

ocTrj.x=linspace(0,1,n);
ocTrj.y=[x0;initialcostate];
ocTrj.y=ocTrj.y(:,ones(1,n));
ocTrj.arcposition=[1;n];
ocTrj.arcinterval=[0 0];
ocTrj.arcarg=arcarg;
ocTrj.x0=0;
ocTrj.timehorizon=0;
ocTrj.modelparameter=parametervalue(ocObj);
ocTrj.modelname=modelname(ocObj);
ocTrj=octrajectory(ocTrj);
