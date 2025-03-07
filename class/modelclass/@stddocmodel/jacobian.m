function varargout=jacobian(docObj,solObj,varargin)
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


if isempty(docObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(docObj);
if ismapprimitive(solObj) || isdoctrajectory(solObj) || isdocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=[initialtime(solObj) independentvar(solObj)];
    depvar=[initialstate(solObj) dependentvar(solObj)];
    arcpos=arcposition(solObj);
    %arcpos(2,end)=arcpos(2,end)-1;
    arcn=arcnum(solObj);
end

if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values 
    varargout{1}=feval(docObj,'SymbolicCanonicalSystemMapJacobian',arcarg);
else

    for ii=1:arcn
        equationnum=canonicalmapdimension(docObj,arcarg(ii));
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        J=zeros(equationnum,2*equationnum,length(arcp)-1);
        for jj=1:length(arcp)-1
            J(:,:,jj)=feval(docObj,'CanonicalSystemMapJacobian',indepvar(arcp(jj+1)),depvar(:,arcp(jj:jj+1)),par,arcarg(ii));
        end
        varargout{1}(1:equationnum,1:2*equationnum,arcp(1:end-1))=J;
    end
end
