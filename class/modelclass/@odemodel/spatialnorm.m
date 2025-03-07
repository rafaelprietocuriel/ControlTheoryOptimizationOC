function varargout=spatialnorm(odeObj,solObj,varargin)
%
% NORM returns the state values/variables.
% 
% X=NORM(OCOBJ) OCOBJ is a stdocmodel class. X is a cell array of strings
% consisting of the state variable names. 
% 
% X=NORM(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for NORM(OCOBJ). Otherwise the state values are
% returned. If SOLOBJ is an octrajectory consisting of multiple arcs X is a
% cell array of matrices, with the state values for each arc separately.
%
% X=NORM(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% NORM(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the state values of all arcs are returned in one
% matrix.
connectflag=[];
if isempty(odeObj)
    return
end
if nargin==1
    solObj=[];
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
try
    N=spacedata(odeObj);
catch
    N=1;
end
if isstruct(solObj)
    try
        depvar=solObj.y;
        arcpos=solObj.arcposition;
        arcn=size(arcpos,2)-1;
    catch
        ocmaterror('If second input argument is a structure the fieldnames ''y'' and ''arcposition'' have to exist!')
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(odeObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isoccurve(solObj)
    arcn=1;
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcpos=[1 size(depvar,2)]';
    end
end
Np1=N+1;

% weights to distinguish between spatially symmetric points
if mod(Np1,2)==0
    wl=ones(1,Np1/2);
    wl=wl/sum(wl);
    w=[wl 1*wl];
else
    wl=ones(1,N/2);
    wl=wl/sum(wl);
    w=[wl 1 1*wl];
end
w=w(:)/sum(w);

n=statenum(odeObj)/Np1;
xcoord=statecoord(odeObj);
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    for jj=1:n
        %x=w(:,ones(1,length(arcp))).*depvar(xcoord((jj-1)*Np1+1:jj*Np1),arcp);
        %l=w(:,ones(1,length(arcp))).*depvar(lcoord((jj-1)*Np1+1:jj*Np1),arcp);
        x=depvar(xcoord((jj-1)*Np1+1:jj*Np1),arcp);
        nrmx=sum(((x(1:N,:)+x(2:Np1,:))/2).^2)/(N+1);
        if connectflag
            varargout{1}(jj,arcp)=sqrt(nrmx);
        else
            varargout{ii}(jj,:)=sqrt(nrmx);
        end
    end
end
