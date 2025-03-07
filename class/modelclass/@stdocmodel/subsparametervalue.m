function s=subsparametervalue(ocObj,s,varargin)
% SUBSPARAMETERVALUE substitutes the model parameters
%
% S=SUBSPARAMETERVALUE(OCOBJ,S) substitutes the parameter values of the
% optimal control model OCOBJ for the expression S, where S can be a
% character or symbolic variable.
%
% S=SUBSPARAMETERVALUE(OCOBJ,S,IDX) substitutes the IDX'th parameter
% values, where IDX can be the position of the parameter value or the
% variable name of the parameter.
%
% S=SUBSPARAMETERVALUE(OCOBJ,S,...,'EXCLUDE',IDX) the IDX'th parameter
% values are excluded from substitution.

if isempty(ocObj)
    return
end
symkernel=getsymkernel();

excludeidx=[];
excludeflag=find(strcmpi(varargin,'exclude'));
if ~isempty(excludeflag)
    excludeidx=varargin{excludeflag+1};
    excludeidx=parameterindex(ocObj,excludeidx);
end

[parval,parvar]=parametervalue(ocObj,varargin{:});
if ~isempty(excludeidx)
    parval(excludeidx)=[];
    parvar(excludeidx)=[];
end
if ~isempty(symkernel)
    if ischar(s)
        s=mystr2sym(s);
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