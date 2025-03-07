function varargout=line(hocTrj,varargin)
%
% LINE basic plot command of an octrajectory.
%
% LINE(OCTRJ,XCOORD,YCOORD) OCTRJ is an instance of the class octrajecotry
% and XCOORD, YCOORD, respectively coordinates of the dependent variables
% of OCTRJ. This commands calls the native line command of MATLAB with:
%   LINE(OCTRJ.Y(XCOORD,:),OCTRJ.Y(YCOORD,:))
% For a detailed description of the LINE command see the MATLAB help.
%
% LINE(OCTRJ,XCOORD,YCOORD,ZCOORD) creates a 3D plot using a third
% coordinate ZCOORD. 
%
% LINE(OCTRJ,XCOORD,YCOORD,(ZCOORD),'PROPERTYNAME',PROPERTYVALUE,...)
% additionally usual line properties can be provided.
%
% LINE(OCTRJ,'XDATA',XDATA,'XCOORD',XCOORD,'YDATA',YDATA,'YCOORD',YCOORD,'Z
% DATA',ZDATA,'ZCOORD',ZCOORD,'PROPERTYNAME',PROPERTYVALUE,...) is the full
% version of the LINE command, with paired arguments. 'XDATA', 'YDATA' and
% optionally 'ZDATA' denote the name of the depicted variables. Possible
% values are:
%   'dependentvar' ... dependent variables of OCTRJ, usually state and
%                      costates
%   'independentvar' ... the independent variable of OCTRJ, usually the
%                      time
% If model specific variables are depicted, e.g., state, control,
% Hamiltonian, etc. the OCOBJ instance of the model has to be provided as
% well, by the paired argument 'OCMODEL',OCOBJ. Possible values for 'XDATA',
% 'YDATA' and optionally 'ZDATA' are: 
%   'independentvar'
%   'time'
%   'dependentvar'
%   'state'
%   'costate',
%   'control'
%   'lagrangemultiplier'
%   'canonicalsystem'
%   'hamiltonian'
%   'userfunction'
%
% LINE(OCTRJ,...,'CONNECT',VAL,...) for VAL=0 (default) each arc is plotted
% separately, for VAL=1 all arcs are combined.
%
% H=LINE(...) returns the line handles H  for the plotted arc(s).

% transform into multiarc path
hocTrj=transform2multiarc(hocTrj);
if nargin<=2 
    ocmaterror('Not enough input arguments.')
end
if isnumeric(varargin{1})
    if isnumeric(varargin{2})
        if nargin>=4 && isnumeric(varargin{3})
            [varargout{1:nargout}]=line(hocTrj,'xcoord',varargin{1},'ycoord',varargin{2},'zcoord',varargin{3},'xdata','dependentvar','ydata','dependentvar','zdata','dependentvar',varargin{4:nargin-1});
            return
        else
            [varargout{1:nargout}]=line(hocTrj,'xcoord',varargin{1},'ycoord',varargin{2},'xdata','dependentvar','ydata','dependentvar',varargin{3:nargin-1});
            return
        end
    else
        ocmaterror('Not enough input arguments.')
    end
end
if mod(nargin-1,2)
    ocmaterror('Line has to be called in name/property value pairs.')
end
modelidx=find(strcmpi(varargin,'ocmodel'));
connectidx=find(strcmpi(varargin,'connect'));
xcoordidx=find(strcmpi(varargin,'xcoord'));
ycoordidx=find(strcmpi(varargin,'ycoord'));
zcoordidx=find(strcmpi(varargin,'zcoord'));
if ~isempty(modelidx)
    ocObj=varargin{modelidx+1};
else
    ocObj=stdocmodel();
end
if ~isempty(connectidx)
    connectflag=varargin{connectidx+1};
else
    connectflag=0;
end
if ~isempty(xcoordidx)
    xcoord=varargin{xcoordidx+1};
else
    xcoord=[];
end
if ~isempty(ycoordidx)
    ycoord=varargin{ycoordidx+1};
else
    ycoord=[];
end
if ~isempty(zcoordidx)
    zcoord=varargin{zcoordidx+1};
else
    zcoord=[];
end
varargin([xcoordidx ycoordidx zcoordidx xcoordidx+1 ycoordidx+1 zcoordidx+1 modelidx modelidx+1 connectidx connectidx+1])=[];
if ~connectflag
    arcn=arcnum(hocTrj);
else
    arcn=1;
end

xdataidx=find(strcmpi(varargin,'xdata'));
ydataidx=find(strcmpi(varargin,'ydata'));
zdataidx=find(strcmpi(varargin,'zdata'));
if isempty(xdataidx)
    ocmaterror('')
end
if isempty(ydataidx)
    ocmaterror('')
end
[xdata{1:arcn}]=propvar2propval(varargin{xdataidx+1});
[ydata{1:arcn}]=propvar2propval(varargin{ydataidx+1});
if ~isempty(zdataidx)
    [zdata{1:arcn}]=propvar2propval(varargin{zdataidx+1});
else
    zdata=[];
end
varargin([xdataidx ydataidx zdataidx xdataidx+1 ydataidx+1 zdataidx+1])=[];

if isempty(xcoord)
    xcoord=size(xdata{1},1);
end
if isempty(ycoord)
    ycoord=size(ydata{1},1);
end
if ~isempty(zdataidx) && isempty(zcoord)
    zcoord=size(zdata{1},1);
end

if ~isempty(zcoord)
    coordnum=[numel(xcoord) numel(ycoord) numel(zcoord)];
else
    coordnum=[numel(xcoord) numel(ycoord)];
end
maxcoordnum=max(coordnum);
mincoordnum=min(coordnum);
%
if maxcoordnum~=mincoordnum
    if mincoordnum>1 || numel(unique(coordnum))>2
        ocmaterror('Vectors must be the same lengths.')
    end
    for ii=1:numel(coordnum)
        if coordnum(ii)==1
            if ii==1
                xcoord=xcoord(1,ones(1,maxcoordnum));
            elseif ii==2
                ycoord=ycoord(1,ones(1,maxcoordnum));
            else
                zcoord=zcoord(1,ones(1,maxcoordnum));
            end
        end
    end
end

counter=0;
h=zeros(arcn*maxcoordnum,1);
arcarg=arcargument(hocTrj);
for ii=1:arcn
    for jj=1:maxcoordnum
        counter=counter+1;
        if ~isempty(zdata)
            lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:),'Zdata',zdata{ii}(zcoord(jj),:)};
        else
            lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:)};
        end
        h(counter)=line(lineargument{:},varargin{:});
        if modelidx
            set(h(counter),'Tag',num2str(arcarg(ii)))
        end
    end
end
if nargout==1
    varargout{1}=h;
end
    function varargout=propvar2propval(propval)
        
        switch lower(propval)
            case {'state'}
                [varargout{1:nargout}]=state(ocObj,hocTrj,connectflag);
            case {'dependentvar','y'}
                y=dependentvar(hocTrj);
                if ~connectflag
                    numberarc=arcnum(hocTrj);
                    arcpos=arcposition(hocTrj);
                    for arc=1:numberarc
                        varargout{arc}=y(:,arcpos(1,arc):arcpos(2,arc));
                    end
                else
                    varargout{1}=y;
                end
            case {'time'}
                [varargout{1:nargout}]=time(ocObj,hocTrj,connectflag);
            case {'impulsetime'}
                [varargout{1:nargout}]=impulsetime(hocTrj);
            case {'independentvar','x'}
                x=independentvar(hocTrj);
                if ~connectflag
                    numberarc=arcnum(hocTrj);
                    arcpos=arcposition(hocTrj);
                    for arc=1:numberarc
                        varargout{arc}=x(arcpos(1,arc):arcpos(2,arc));
                    end
                else
                    varargout{1}=x;
                end
            case {'costate'}
                [varargout{1:nargout}]=costate(ocObj,hocTrj,connectflag);
            case {'constraint'}
                [varargout{1:nargout}]=constraint(ocObj,hocTrj,connectflag);
            case {'impulsecontrol'}
                [varargout{1:nargout}]=impulsecontrol(ocObj,hocTrj);
            case {'control'}
                [varargout{1:nargout}]=control(ocObj,hocTrj,connectflag);
            case {'hamiltonian'}
                [varargout{1:nargout}]=hamiltonian(ocObj,hocTrj,connectflag);
            case {'lagrangemultiplier'}
                [varargout{1:nargout}]=lagrangemultiplier(ocObj,hocTrj,connectflag);
            case {'userfunction'}
                [varargout{1:nargout}]=userfunction(ocObj,hocTrj,connectflag);
            case {'objectivefunction'}
                [varargout{1:nargout}]=objectivefunction(ocObj,hocTrj,connectflag);
            case {'canonicalsystem'}
                [varargout{1:nargout}]=canonicalsystem(ocObj,hocTrj,connectflag);
            case {'adjointsystem'}
                [varargout{1:nargout}]=adjointsystem(ocObj,hocTrj);
            case {'spatialnorm'}
                [varargout{1:nargout}]=spatialnorm(ocObj,hocTrj);
            otherwise
                [varargout{1:nargout}]=propval;
        end
    end
end