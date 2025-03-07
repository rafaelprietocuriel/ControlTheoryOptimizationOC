function J=linearize(dynPrim,ocObj,lineartype,varargin)
%
J=[];
dynprimitiveflag=0;
if isempty(ocObj) || isempty(dynPrim)
    return
end
if nargin>=4
    dynprimitiveflag=varargin{1};
end
if isequilibrium(dynPrim) || isfixpoint(dynPrim)
    switch lineartype
        case 'dependentvar'
            J=jacobian(ocObj,dynPrim,[],dynprimitiveflag,varargin{2:end});
    end
else
    ocmatmsg('Dynprimitive class ''%s'' is unknown.\n',dynprimitiveclass(dynPrim))
end

