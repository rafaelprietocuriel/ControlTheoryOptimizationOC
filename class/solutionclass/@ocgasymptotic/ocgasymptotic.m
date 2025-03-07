function ocgAsym=ocgasymptotic(varargin)
%
% OCGASYMPTOTIC ocasymptotic constructor
%
idx=find(cellfun('isempty',varargin));
if ~isempty(idx)
    varargin(idx)=[];
    ocgAsym=ocgasymptotic(varargin{:});
    return
end

switch nargin
    case 0
        ocgAsym.limitset=limitset(gdynprimitive());
        ocgAsym=class(ocgAsym,'ocgasymptotic',ocgtrajectory(),gdynprimitive());

    case 1
        if isa(varargin{1},'ocgasymptotic')
            ocgAsym=varargin{1};
        else
            error('Wrong input argument.')
        end
    case 2
        idx=find(cellfun('isclass',varargin,'ocgtrajectory'));
        if isempty(idx)
            error('No ''ocgtrajectory''.')
        end
        ocgTrj=varargin{idx};
        varargin(idx)=[];
        idx=find(cellfun('isclass',varargin,'gdynprimitive'));
        if isempty(idx)
            error('No ''gdynprimitive''.')
        end
        ocgAsym.limitset=varargin{idx};
        ocgAsym=class(ocgAsym,'ocgasymptotic',ocgTrj,varargin{idx});
end