function s=subsparametervalue(odeObj,s,varargin)
% SUBSPARAMETERVALUE substitutes the model parameters
%
% S=SUBSPARAMETERVALUE(ODEOBJ,S) substitutes the parameter values of the
% optimal control model ODEOBJ for the expression S, where S can be a
% character or symbolic variable.

if isempty(odeObj)
    return
end
symkernel=getsymkernel();

excludeidx=[];
excludeflag=find(strcmpi(varargin,'exclude'));
if ~isempty(excludeflag)
    excludeidx=varargin{excludeflag+1};
    excludeidx=parameterindex(odeObj,excludeidx);
end

[parval,parvar]=parametervalue(odeObj,varargin{:});
if ~isempty(excludeidx)
    parval(excludeidx)=[];
    parvar(excludeidx)=[];
end
if ~isempty(symkernel)
    if ischar(s)
        s=sym(s);
    end
    if ~isa(s,'sym')
        ocmaterror('Second argument is not a symbolic expression.')
    end
    switch symkernel
        case 'mupad'
            parval=num2cell(parval);
    end
    s=subs(s,parvar,parval);
else
    if isa(s,'sym')
        s=char(s);
    end

    ocmatmsg('Symbolic toolbox not activated. Substitution is done as string replacement.\n')
    for ii=1:numel(parval)
        s=regexprep(s,['\<' parvar{ii} '\>'],num2str(parval(ii),10));
    end
end