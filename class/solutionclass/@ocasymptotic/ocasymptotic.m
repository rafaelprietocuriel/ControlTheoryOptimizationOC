function ocAsym=ocasymptotic(varargin)
%
% OCASYMPTOTIC ocasymptotic constructor
%
idx=find(cellfun('isempty',varargin));
if ~isempty(idx)
    varargin(idx)=[];
    ocAsym=ocasymptotic(varargin{:});
    return
end

switch nargin
    case 0
        ocAsym.limitset=limitset(dynprimitive());
        ocAsym=class(ocAsym,'ocasymptotic',octrajectory(),dynprimitive());
    case 1
        if isocasymptotic(varargin{1})
            ocAsym=varargin{1};
        else
            error('Input argument is not an ''ocasymptotic''.')
        end
    case 2
        idx=find(cellfun('isclass',varargin,'octrajectory'));
        if isempty(idx)
            error('No ''octrajectory''.')
        end
        ocTrj=varargin{idx};
        varargin(idx)=[];
        idx=find(cellfun('isclass',varargin,'dynprimitive'));
        if isempty(idx)
            error('No ''dynprimitive''.')
        end
        ocAsym.limitset=varargin{idx};
        ocAsym=class(ocAsym,'ocasymptotic',ocTrj,varargin{idx});
end