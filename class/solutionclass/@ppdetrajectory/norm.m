function out=norm(ppdeTrj,varargin)

if ~ischar(varargin{1})
    ocmaterror('Second argument is not a string.')
end
switch lower(varargin{1})
    case 'l2state'
        out=norm(ppdePrim.ppdetrajectory,varargin{:});
    case 'l2costate'
end
