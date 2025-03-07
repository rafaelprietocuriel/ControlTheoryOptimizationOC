function out=time(ocgAsym,connectflag)
if nargin==1
    connectflag=1;
end
ocObj=loadmodel(ocgAsym);

out=time(ocObj,ocgAsym.ocgtrajectory,connectflag);
