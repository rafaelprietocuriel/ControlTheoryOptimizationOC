function out=time(ocgTrj,connectflag)
if nargin==1
    connectflag=1;
end
ocObj=loadmodel(ocgTrj);

out=time(ocObj,ocgTrj.octrajectory,connectflag);
