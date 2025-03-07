function s=subsparametervalue(ocObj,s)
% SUBSPARAMETERVALUE substitutes the model parameters
%
% S=SUBSPARAMETERVALUE(OCOBJ,S) substitutes the parameter values of the
% optimal control model OCOBJ for the expression S, where S can be a
% character or symbolic variable.

if isempty(ocObj)
    return
end
symkernel=getsymkernel();

[parval,parvar]=parametervalue(ocObj);
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