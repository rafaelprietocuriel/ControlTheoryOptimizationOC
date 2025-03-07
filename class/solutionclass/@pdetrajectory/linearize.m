function J=linearize(pdeTrj,ppdeObj,lineartype,varargin)
%
J=[];
if isempty(ppdeObj) || isempty(pdeTrj)
    return
end
if isequilibrium(pdeTrj)
    switch lineartype
        case 'dependentvar'
            J=jacobian(ppdeObj,pdeTrj,[],varargin{:});
    end
else
    ocmatmsg('Ppdeprimitive class ''%s'' is unknown.\n',ppdeprimitiveclass(pdeTrj))
end

