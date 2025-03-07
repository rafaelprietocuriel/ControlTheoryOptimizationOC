function out=solverinfo(ocTrj,varargin)

if nargin>1
    try
        out=ocTrj.solverinfo.(varargin{1});
        return
    end
end
out=ocTrj.solverinfo;
