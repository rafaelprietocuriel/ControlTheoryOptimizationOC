function J=linearize(ppdePrim,ppdeObj,lineartype,varargin)
%
J=[];
dynprimitiveflag=0;
if isempty(ppdeObj) || isempty(ppdePrim)
    return
end
if nargin>=4
    dynprimitiveflag=varargin{1};
end
if isequilibrium(ppdePrim)
    switch lineartype
        case 'dependentvar'
            J=jacobian(ppdeObj,ppdePrim,[],dynprimitiveflag,varargin{2:end});
    end
else
    ocmatmsg('Ppdeprimitive class ''%s'' is unknown.\n',ppdeprimitiveclass(ppdePrim))
end

