function varargout=jacobian(dgObj,varargin)
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


solstruct=[];
arcarg=[];
dynprimitiveflag=0;
jacobianfile='';
if isempty(dgObj)
    return
end
if nargin==1
    solstruct=[];
end
if nargin>=2
    solstruct=varargin{1};
end
if nargin>=3
    arcarg=varargin{2};
end
if nargin>=4
    dynprimitiveflag=varargin{3};
end
jacobianfileidx=find(strcmpi(varargin,'jacobianfile'));
if ~isempty(jacobianfileidx)
    jacobianfile=varargin{jacobianfileidx+1};
end
if isempty(jacobianfile)
    jacobianfile='CanonicalSystemJacobian';
end
if isstruct(solstruct)
    try
        arcarg=solstruct.arcarg;
        indepvar=solstruct.x;
        depvar=solstruct.y;
        arcn=1;
        arcpos=[1;1];
    catch
        ocmaterror('If the second input argument is a structure the fields ''y'', ''x'' and ''arcarg'' have to exist!')
    end
elseif isdynprimitive(solstruct) || isoctrajectory(solstruct) || isocasymptotic(solstruct) || isoccurve(solstruct)
    arcarg=arcargument(solstruct);
    indepvar=time(dgObj,solstruct,1);
    depvar=dependentvar(solstruct);
    arcpos=arcposition(solstruct);
    arcn=arcnum(solstruct);
end
if isempty(arcarg)
    arcarg=0;
end

if isempty(solstruct)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values 
    varargout{1}=feval(dgObj,'SymbolicCanonicalSystemJacobian',arcarg);
else
    % return optimal control value evaluated at 'solstruct'
    %     if nargin>=4
    %         indepvar=varargin{1};
    %     end
    if dynprimitiveflag
        try
            %varargout{1}=feval(dgObj,'CanonicalSystemJacobian',indepvar,depvar,parametervalue(dgObj),arcarg);
            if strcmp(jacobianfile,'CanonicalSystemJacobian')
                varargout{1}=feval(dgObj,jacobianfile,indepvar,depvar,parametervalue(dgObj),arcarg);
            else
                varargout{1}=feval(dgObj,jacobianfile,depvar,parametervalue(dgObj),arcarg);
            end
        catch
            varargout{1}=feval(dgObj,'CanonicalSystemJacobian',indepvar,depvar,parametervalue(dgObj),arcarg);
        end
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
            equationnum=canonicalsystemdimension(dgObj);
            J=zeros(equationnum*equationnum,size(arcp,2));
            for jj=1:size(arcp,2)
                tmp=feval(dgObj,'CanonicalSystemJacobian',indepvar(arcp(jj)),depvar(:,arcp(jj)),parametervalue(dgObj),arcarg(ii));
                if size(arcp,2)>1
                    J(:,jj)=tmp(:);
                else
                    J=tmp;
                end
            end
            if connectflag
                varargout{1}(1:equationnum*equationnum,arcp)=J;
            else
                varargout{ii}=J;
            end
        end
    end
end
