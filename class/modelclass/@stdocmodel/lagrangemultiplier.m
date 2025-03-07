function varargout=lagrangemultiplier(ocObj,solObj,varargin)
%
%
% LAGRANGEMULTIPLIER returns the values of the Lagarangian multipliers.
% 
% LM=LAGRANGEMULTIPLIER(OCOBJ) OCOBJ is a stdocmodel class. LM is a cell
% array of strings consisting of the terms of the Lagarangian multiplier
% derived from the maximum Hamiltonian condition (explicit case) or the
% variable names of Lagarangian multipliers (implicit case).  
% 
% LM=LAGRANGEMULTIPLIER(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for LAGRANGEMULTIPLIER(OCOBJ). Otherwise the values
% of the Lagarangian multipliers are returned. If SOLOBJ is an octrajectory
% consisting of multiple arcs LM is a cell array of matrices, with the
% values of the Lagarangian multipliers for each arc separately. 
%
% LM=LAGRANGEMULTIPLIER(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% LAGRANGEMULTIPLIER(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are
% returned "connected", i.e., values of the Lagarangian multipliers of all
% arcs are returned in one matrix.

if isempty(ocObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
        arcpos=solObj.arcposition;
        arcn=numel(arcarg);
    catch
        ocmaterror('If the second input argument is a structure the fields ''dependentvar'', ''independentvar'' and ''arcarg'' have to exist!')
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif  isoccurve(solObj)
    arcarg=arcargument(solObj);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcarg=unique(arcarg);
        arcpos=[1 size(depvar,2)];
        arcn=1;
    else
        arcpos=arcposition(solObj);
        try
            arcn=arcnum(solObj);
        catch
            arcn=1;
        end
    end
    indepvar=zeros(1,size(depvar,2));
end
if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values
    repflag=[];
    arcarg=[];
    if nargin>=3
        arcarg=varargin{1};
    end
    if isempty(arcarg)
        varargout{1}=lagrangevariable(ocObj);
    else
        [varargout{1:nargout}]=feval(ocObj,'SymbolicLagrangeMultiplier',arcarg);
    end
else
    % return value of the canonical system evaluated at 'depvar'
    connectflag=[];
    if nargin>=3
        connectflag=varargin{1};
    end
    if isempty(connectflag)
        connectflag=0;
    end
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        [lmmc,lmsc]=feval(ocObj,'LagrangeMultiplier',indepvar(arcp),depvar(:,arcp),par,arcarg(ii));
        if connectflag
            varargout{1}(:,arcp)=[lmmc;lmsc];
        else
            varargout{ii}=[lmmc;lmsc];
        end
    end
end
