function varargout=userfunction(docObj,solObj,varargin)
%
% CONTROL returns the control values.
% 
% U=CONTROL(OCOBJ) OCOBJ is a stdocmodel class. U is a cell array of
% strings consisting of the terms derived from the maximum Hamiltonian
% condition (explicit control) or control variable names (implicit
% control).  
% 
% U=CONTROL(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for CONTROL(OCOBJ). Otherwise the control values
% are returned. If SOLOBJ is an octrajectory consisting of multiple arcs U
% is a cell array of matrices, with the control values for each arc
% separately. 
%
% U=CONTROL(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% CONTROL(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the control values of all arcs are returned in one
% matrix.
if isempty(docObj)
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(docObj);
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
    catch
        ocmaterror('If the second input argument is a structure the fields ''dependentvar'', ''independentvar'' and ''arcarg'' have to exist!')
    end
elseif ismapprimitive(solObj) || isdoctrajectory(solObj) || isdocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=[initialtime(solObj) independentvar(solObj)];
    depvar=[initialstate(solObj) dependentvar(solObj)];
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
end

if isempty(solObj)
    try
        varargout{1}=feval(docObj,'UserFunction',[],[],par,[]);
    catch
        varargout{1}=[];
        return
    end
else
    % return value of the canonical system evaluated at 'depvar'
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        out=feval(docObj,'UserFunction',indepvar(arcp(2:end)),depvar(:,arcp),par,arcarg(ii));
        varargout{1}(1:ctrlnum,arcp(1:end-1))=out;
    end
end
