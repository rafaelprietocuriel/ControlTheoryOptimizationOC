function [ocgTrjNew,tnew,errval]=comptimegrid(ocgTrj,ocgTrjFine,varargin)
% TIME returns the time argument of the ocgradtrajectory
mindist=[];
gamma=[];
opt=[];
if isempty(modelname(ocgTrj))
    errval=[];
    return
end
mindistidx=find(strcmpi(varargin,'mindist'));
gammaidx=find(strcmpi(varargin,'gamma'));
optionidx=find(strcmpi(varargin,'option'));
if ~isempty(mindistidx)
    mindist=varargin{mindistidx+1};
end
if ~isempty(gammaidx)
    gamma=varargin{gammaidx+1};
end
if ~isempty(optionidx)
    opt=varargin{optionidx+1};
end
if isempty(mindist)
    mindist=1e-1;
end
if isempty(gamma)
    gamma=1e3;
end
if isempty(opt)
    opt=setocoptions('EQ','LargeScale','off');
end

ocObj=stdocmodel(modelname(ocgTrj));
par=modelparameter(ocgTrj);
t=time(ocgTrj);
n=length(t);

h=diff(t);
tf=time(ocgTrjFine);
funch=modelspecificfunc(ocObj,'4GradContinuation');
funch=funch();
statedynamics=funch{1}{1};
costatedynamics=funch{2}{1};

dXdt=[statedynamics(t,state(ocgTrj),control(ocgTrj),par); ...
    costatedynamics(t,costate(ocgTrj),[state(ocgTrj);control(ocgTrj)],par)];

dXdtf=[statedynamics(tf,state(ocgTrjFine),control(ocgTrjFine),par); ...
    costatedynamics(tf,costate(ocgTrjFine),[state(ocgTrjFine);control(ocgTrjFine)],par)];

dXdtmean=repmat(h,size(dXdt,1),[])/2.*(dXdt(:,2:end)+dXdt(:,1:end-1));


errval=dXdtmean-dXdtf(:,2:2:end-1);
normerrval=sqrt(sum(errval.^2));
hnew=quadprog(1,gamma*normerrval.',[],[],ones(1,n-1),t(end),repmat(mindist,n-1,1),repmat(t(end),n-1,1),[],opt.EQ);
tnew=cumsum([t(1) hnew.']);
tnew(end)=t(end);
ocgTrjNew=deval(ocgTrj,tnew);
plot(tf(2:2:end-1),normerrval,tf(2:2:end-1),hnew)