function varargout=line(ocgTrj,varargin)
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
if nargin<=2 
    ocmaterror('Not enough input arguments.')
end
if isnumeric(varargin{1})
    if isnumeric(varargin{2})
        if nargin>=4 && isnumeric(varargin{3})
            [varargout{1:nargout}]=line(ocgTrj,'xcoord',varargin{1},'ycoord',varargin{2},'zcoord',varargin{3},'xdata','state','ydata','state','zdata','state',varargin{4:nargin-1});
            return
        else
            [varargout{1:nargout}]=line(ocgTrj,'xcoord',varargin{1},'ycoord',varargin{2},'xdata','state','ydata','costate',varargin{3:nargin-1});
            return
        end
    else
        ocmaterror('Not enough input arguments.')
    end
end
if mod(nargin-1,2)
    ocmaterror('Line has to be called in name/property value pairs.')
end
xcoordidx=find(strcmpi(varargin,'xcoord'));
ycoordidx=find(strcmpi(varargin,'ycoord'));
zcoordidx=find(strcmpi(varargin,'zcoord'));
modelidx=find(strcmpi(varargin,'ocmodel'));
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
varargin([xcoordidx ycoordidx zcoordidx xcoordidx+1 ycoordidx+1 zcoordidx+1 modelidx modelidx+1])=[];
arcn=1;

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
for ii=1:arcn
    for jj=1:maxcoordnum
        counter=counter+1;
        if ~isempty(zdata)
            lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:),'Zdata',zdata{ii}(zcoord(jj),:)};
        else
            lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:)};
        end
        h(counter)=line(lineargument{:},varargin{:});
    end
end
if nargout==1
    varargout{1}=h;
end
    function varargout=propvar2propval(propval)
        
        switch lower(propval)
            case {'state'}
                [varargout{1:nargout}]=state(ocgTrj);
            case {'y'}
                y=dependentvar(ocgTrj);
                if ~connectflag
                    numberarc=1;
                    arcpos=arcposition(ocgTrj);
                    for arc=1:numberarc
                        varargout{arc}=y(:,arcpos(1,arc):arcpos(2,arc));
                    end
                else
                    varargout{1}=y;
                end
            case {'time'}
                [varargout{1:nargout}]=time(ocgTrj);
            case {'t'}
                x=independentvar(ocgTrj);
                if ~connectflag
                    numberarc=1;
                    arcpos=arcposition(ocgTrj);
                    for arc=1:numberarc
                        varargout{arc}=x(arcpos(1,arc):arcpos(2,arc));
                    end
                else
                    varargout{1}=x;
                end
            case {'modelparameter'}
                par=modelparameter(ocgTrj);
                par=par(:);
                arcpos=arcposition(ocgTrj);
                if ~connectflag
                    numberarc=arcnum(ocgTrj);
                    for arc=1:numberarc
                        varargout{arc}=par(:,ones(1,arcpos(2,arc)-arcpos(1,arc)+1));
                    end
                else
                    varargout{1}=par(:,ones(1,arcpos(2,end)));
                end
            case {'costate'}
                [varargout{1:nargout}]=costate(ocgTrj);
            case {'constraint'}
                [varargout{1:nargout}]=constraint(ocgTrj);
            case {'control'}
                [varargout{1:nargout}]=control(ocgTrj);
            case {'hamiltonian'}
                [varargout{1:nargout}]=hamiltonian(ocgTrj);
            case {'userfunction'}
                [varargout{1:nargout}]=userfunction(ocgTrj);
            case {'objectivefunction'}
                [varargout{1:nargout}]=objectivefunction(ocgTrj);
            case {'canonicalsystem'}
                [varargout{1:nargout}]=canonicalsystem(ocgTrj);
            case {'adjointsystem'}
                [varargout{1:nargout}]=adjointsystem(ocgTrj);
            otherwise
                [varargout{1:nargout}]=propval;
        end
    end
end