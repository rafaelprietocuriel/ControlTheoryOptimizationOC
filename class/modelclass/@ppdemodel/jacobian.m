function varargout=jacobian(ppdeObj,varargin)
%
% JACOBIAN returns the symbolic Jacobian or evaluated at a solution object.
%
% J=JACOBIAN(OCOBJ) OCOBJ is a stdocmodel class. J is a string matrix of
% the formally written Jacobian for the canonical system.
%
% J=JACOBIAN(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for JACOBIAN(OCOBJ). Otherwise the Jacobian
% evaluated at the solution object is returned. If SOLOBJ is an
% octrajectory consisting of multiple arcs J is a cell array of
% 3D-matrices, with the Jacobian values for each arc separately.
%
% J=JACOBIAN(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% JACOBIAN(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the Jacobian values of all arcs are returned in one
% 3D-matrix.


J=[];
solObj=[];
arcarg=[];
if isempty(ppdeObj)
    return
end
if nargin==1
    solObj=[];
end
if nargin>=2
    solObj=varargin{1};
end
if nargin>=3
    arcarg=varargin{2};
end
if nargin>=4
    femdat=varargin{3};
end
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.x;
        depvar=solObj.y;
        arcn=1;
        arcpos=[1;1];
    catch
        ocmaterror('If the second input argument is a structure the fields ''y'', ''x'' and ''arcarg'' have to exist!')
    end
elseif ispdeprimitive(solObj) || ispdetrajectory(solObj) || ispdeasymptotic(solObj)
    femdat=femdata(solObj);
    arcarg=arcargument(solObj);
    indepvar=time(ppdeObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
end
if isempty(arcarg)
    arcarg=0;
end

if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values
    varargout{1}=[];
else
    % return optimal control value evaluated at 'solObj'
    %     if nargin>=4
    %         indepvar=varargin{1};
    %     end
    if length(indepvar)==1
        varargout{1}=feval(ppdeObj,'CanonicalSystemJacobian',indepvar,depvar,parametervalue(ppdeObj),arcarg,femdat);
    else
        connectflag=[];
        if nargin>=4
            connectflag=varargin{4};
        end
        if isempty(connectflag)
            connectflag=0;
        end
        for ii=1:arcn
            arcp=arcpos(1,ii):arcpos(2,ii);
            %equationnum=canonicalsystemdimension(ppdeObj,arcarg(ii));
            %J=zeros(equationnum*equationnum,size(arcp,2));
            for jj=1:size(arcp,2)
                tmp=feval(ppdeObj,'CanonicalSystemJacobian',indepvar(arcp(jj)),depvar(:,arcp(jj)),parametervalue(ppdeObj),arcarg(ii),femdat);
                if size(arcp,2)>1
                    J(:,jj)=tmp(:);
                else
                    J=tmp;
                end
            end
            if connectflag
                varargout{1}(:,:,arcp)=J;
            else
                varargout{ii}=J;
            end
        end
    end
end
