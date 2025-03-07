function varargout=propvar2propval(ocObj,ocTrj,propval,varargin)
if nargin==3
    connectflag=0;
else
    connectflag=varargin{1};
end
switch lower(propval)
    case {'state'}
        [varargout{1:nargout}]=state(ocObj,ocTrj,connectflag);
    case {'dependentvar','y'}
        y=dependentvar(ocTrj);
        if nargout>1
            numberarc=arcnum(ocTrj);
            arcpos=arcposition(ocTrj);
            for arc=1:numberarc
                varargout{arc}=y(:,arcpos(1,arc):arcpos(2,arc));
            end
        else
            varargout{1}=y;
        end
    case {'time'}
        [varargout{1:nargout}]=time(ocObj,ocTrj,connectflag);
    case {'independentvar','x'}
        [varargout{1:nargout}]=independentvar(ocTrj);
    case {'costate'}
        [varargout{1:nargout}]=costate(ocObj,ocTrj,connectflag);
    case {'control'}
        [varargout{1:nargout}]=control(ocObj,ocTrj,connectflag);
    case {'hamiltonian'}
        [varargout{1:nargout}]=hamiltonian(ocObj,ocTrj,connectflag);
    case {'lagrangemultiplier'}
        [varargout{1:nargout}]=lagrangemultiplier(ocObj,ocTrj,connectflag);
    case {'userfunction'}
        [varargout{1:nargout}]=userfunction(ocObj,ocTrj);
    case {'canonicalsystem'}
        [varargout{1:nargout}]=canonicalsystem(ocObj,ocTrj,connectflag);
    case {'adjointsystem'}
        [varargout{1:nargout}]=adjointsystem(ocObj,ocTrj);
    otherwise
        [varargout{1:nargout}]=propval;
end
