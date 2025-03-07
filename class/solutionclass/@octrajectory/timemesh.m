function varargout=timemesh(varargin)
%

zdata=[];
xdata='';
axhdl=[];

if ishandle(varargin{1})
    if ~strcmp(get(varargin{1},'type'),'axes')
        return
    end
    axhdl=varargin{1};
    varargin(1)=[];
end
modelidx=find(strcmpi(varargin,'ocmodel'));
if ~isempty(modelidx)
    ocObj=varargin{modelidx+1};
else
    ocObj=stdocmodel();
end
ocTrj=varargin{1};
varargin(1)=[];
if ~isoctrajectory(ocTrj)
    ocmaterror('First or second argument is not an ''octrajectory''.')
end

if isempty(ocTrj)
    if nargout==1
        varargout{1}=[];
    end
    return
end

if length(varargin)<2
    ocmaterror('Not enough input arguments')
end
xcoord=varargin{1};
if ~isnumeric(xcoord) || any(rem(xcoord,1)) || any(xcoord<0)
    ocmaterror('X coordinate must be a vector of nonnegative integers.')
end
varargin(1)=[];

if ~isempty(varargin)
    if isnumeric(varargin{1})
        zcoord=varargin{1};
        if  any(rem(zcoord,1)) || any(zcoord<0)
            ocmaterror('Z coordinate must be a vector of nonnegative integers.')
        else
            varargin(1)=[];
        end
    end
end

xdataidx=find(strcmpi(varargin,'xdata'));
zdataidx=find(strcmpi(varargin,'zdata'));

remidx=[];
if isempty(xcoord)
    ocmaterror('Vectors must be the same lengths.')
elseif isempty(xcoord) && isempty(y)
    [varargout{1:nargout}]=plot([],[]);
    return
end

if ~isempty(xdata) && ~isempty(xdataidx)
    ocmatmsg('XData is set to independent variable.')
    varargin(xdataidx:xdataidx+1)=[];
    xdataidx=[];
end
if isempty(xdata) && isempty(xdataidx)
    xdata='dependentvar';
elseif ~isempty(xdataidx)
    xdata=varargin{xdataidx+1};
    remidx=[remidx xdataidx:xdataidx+1];
end
if isempty(zdataidx)
    zdata='dependentvar';
else
    zdata=varargin{zdataidx+1};
    remidx=[remidx zdataidx:zdataidx+1];
end
varargin(remidx)=[];
arcn=1;
connectflag=0;
[xdataval{1:arcn}]=propvar2propval(xdata);
[ydataval{1:arcn}]=propvar2propval('time');
[zdataval{1:arcn}]=propvar2propval(zdata);

for ii=1:arcn
    xnum=length(xdataval{ii});
    ynum=length(ydataval{ii});
    X=repmat(xdataval{ii},1,ynum);
    Y=repmat(ydataval{ii},xnum,1);
    Z=zdataval{ii}(zcoord,:);
    h=plot3(X,Y,Z);
end

if nargout==1
    varargout{1}=h;
end
    function varargout=propvar2propval(propval)
        
        switch lower(propval)
            case {'state'}
                [varargout{1:nargout}]=state(ocObj,ocTrj,connectflag);
            case {'extendedstate'}
                [varargout{1:nargout}]=extendedstate(ocObj,ocTrj,connectflag);
            case {'dependentvar','y'}
                y=dependentvar(ocTrj);
                if ~connectflag
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
                x=independentvar(ocTrj);
                if ~connectflag
                    numberarc=arcnum(ocTrj);
                    arcpos=arcposition(ocTrj);
                    for arc=1:numberarc
                        varargout{arc}=x(arcpos(1,arc):arcpos(2,arc));
                    end
                else
                    varargout{1}=x;
                end
            case {'costate'}
                [varargout{1:nargout}]=costate(ocObj,ocTrj,connectflag);
            case {'constraint'}
                [varargout{1:nargout}]=constraint(ocObj,ocTrj,connectflag);
            case {'control'}
                [varargout{1:nargout}]=control(ocObj,ocTrj,connectflag);
            case {'hamiltonian'}
                [varargout{1:nargout}]=hamiltonian(ocObj,ocTrj,connectflag);
            case {'lagrangemultiplier'}
                [varargout{1:nargout}]=lagrangemultiplier(ocObj,ocTrj,connectflag);
            case {'userfunction'}
                [varargout{1:nargout}]=userfunction(ocObj,ocTrj,connectflag);
            case {'canonicalsystem'}
                [varargout{1:nargout}]=canonicalsystem(ocObj,ocTrj,connectflag);
            case {'dynamics'}
                [varargout{1:nargout}]=dynamics(ocObj,ocTrj,connectflag);
            case {'adjointsystem'}
                [varargout{1:nargout}]=adjointsystem(ocObj,ocTrj);
            case 'space'
                L=parametervalue(ocObj,'L');
                [varargout{1:nargout}]=linspace(-L,L,length(xcoord)).';
            otherwise
                [varargout{1:nargout}]=propval;
        end
    end
end
