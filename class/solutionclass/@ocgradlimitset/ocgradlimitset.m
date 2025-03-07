function dynPrim=ocgradlimitset(varargin)
%
%

switch nargin
    case 0
        dynPrim.period=[];
        dynPrim.linearization=[];
        ocgTrj=ocgradtrajectory([]);
        dynPrim=class(dynPrim,'ocgradlimitset',ocgTrj);

    case 1
        if isstruct(varargin{1})
            try
                dynPrim=varargin{1};
                dynPrim=rmfield(dynPrim,'ocgradtrajectory');
                ocgTrj=ocgradtrajectory(varargin{1}.ocgradtrajectory);
                dynPrim=class(dynPrim,'ocgradlimitset',ocgTrj);
            catch
                rethrow(lasterror)
            end
        elseif isocgradlimitset(varargin{1})
            if isempty(varargin{1})
                dynPrim=ocgradlimitset();
            else
                dynPrim=varargin{1};
            end
        end
end
