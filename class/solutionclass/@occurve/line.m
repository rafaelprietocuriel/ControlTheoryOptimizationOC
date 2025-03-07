function varargout=line(ocCuv,varargin)
%
% LINE basic plot command of an octrajectory.
%
% LINE(OCTRJ,'XVAR',XVAR,'XCOORD',XCOORD,'YVAR',YVAR,'YCOORD',YCOORD,'OCMODEL',OCOBJ)
%
% further properties are: 'OCModel'
%           'Connect' ... determines if line is plotted seperately for
%           every arc 

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
    arcn=1;%arcnum(ocCuv);
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
for ii=1:arcn
    for jj=1:maxcoordnum
        counter=counter+1;
        if ~isempty(zdata)
            lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:),'Zdata',zdata{ii}(zcoord(jj),:)};
        else
            lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:)};
        end
        try
            h(counter)=line(lineargument{:},varargin{:});
        catch
            if ~isempty(zdata)
                xdata{ii}=real(xdata{ii});
                ydata{ii}=real(ydata{ii});
                zdata{ii}=real(zdata{ii});
                lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:),'Zdata',zdata{ii}(zcoord(jj),:)};
            else
                xdata{ii}=real(xdata{ii});
                ydata{ii}=real(ydata{ii});
                lineargument={'Xdata',xdata{ii}(xcoord(jj),:),'Ydata',ydata{ii}(ycoord(jj),:)};
            end
            h(counter)=line(lineargument{:},varargin{:});
        end
    end
end
if nargout==1
    varargout{1}=h;
end
    function varargout=propvar2propval(propval)
        switch lower(propval)
            case {'state'}
                [varargout{1:nargout}]=state(ocObj,ocCuv,connectflag);
            case {'constraint'}
                [varargout{1:nargout}]=constraint(ocObj,ocCuv,connectflag);
            case {'dependentvar','y'}
                [varargout{1:nargout}]=dependentvar(ocCuv);
            case {'time'}
                [varargout{1:nargout}]=time(ocObj,ocCuv,connectflag);
            case {'bifpar'}
                [varargout{1:nargout}]=bifurcationparameter(ocCuv);
            case {'independentvar','x'}
                [varargout{1:nargout}]=independentvar(ocCuv);
            case {'costate'}
                [varargout{1:nargout}]=costate(ocObj,ocCuv,connectflag);
            case {'control'}
                [varargout{1:nargout}]=control(ocObj,ocCuv,connectflag);
            case {'hamiltonian'}
                [varargout{1:nargout}]=hamiltonian(ocObj,ocCuv);
            case {'lagrangemultiplier'}
                [varargout{1:nargout}]=lagrangemultiplier(ocObj,ocCuv);
            case {'userfunction'}
                [varargout{1:nargout}]=userfunction(ocObj,ocCuv);
            case {'canonicalsystem'}
                [varargout{1:nargout}]=canonicalsystem(ocObj,ocCuv);
            case {'adjointsystem'}
                [varargout{1:nargout}]=adjointsystem(ocObj,ocCuv);
            case {'spatialnorm'}
                [varargout{1:nargout}]=spatialnorm(ocObj,ocCuv);
            case {'objectivevalue'}
                if isfield(ocCuv.userinfo,'objectivevalue')
                    [varargout{1:nargout}]=ocCuv.userinfo.objectivevalue;
                else
                    [varargout{1:nargout}]=propval;
                end
                
            otherwise
                [varargout{1:nargout}]=propval;
        end
    end
end