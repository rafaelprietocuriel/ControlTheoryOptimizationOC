function ocTrj=horzcat(varargin)
%
% VERTCAT concatenates octrajectories

varargin(cellfun('isempty',varargin)==1)=[];
if numel(varargin)==1
    ocTrj=varargin{1};
    return
end

if numel(varargin)>2
    ocTrj=vertcat(varargin{2:end});
    return
end

asymflag=isocasymptotic(varargin{:});
if sum(asymflag)==2
    ocmaterror('Cannot concatenate two ocasymptotics.')
end
if asymflag(1)
    ocmaterror('Order of octrajectory and ocasymptotic is wrong.')
end
if ~all(isoctrajectory(varargin{:}))
    ocmaterror('Only octrajectories can be concatenated.')
end
y1=dependentvar(varargin{1});
y2=dependentvar(varargin{2});
t1=independentvar(varargin{1});
t2=independentvar(varargin{2});
arcarg1=arcargument(varargin{1});
arcarg2=arcargument(varargin{2});
arcpos1=arcposition(varargin{1});
arcpos2=arcposition(varargin{2});
arcintv1=arcinterval(varargin{1});
arcintv2=arcinterval(varargin{2});

[row1y col1y]=size(y1);
[row2y col2y]=size(y2);
row=max([row1y,row2y]);

trj.y=zeros(row,col1y+col2y);
trj.y(1:row1y,1:col1y)=y1;
trj.y(1:row2y,col1y+(1:col2y))=y2;

trj.x=[t1 t1(end)+t2];
trj.arcarg=[arcarg1 arcarg2];
trj.arcposition=[arcpos1 arcpos1(end)+arcpos2];
trj.arcinterval=[arcintv1 arcintv1(end)+arcintv2(2:end)];
trj.timehorizon=trj.arcinterval(end);
trj.x0=initialtime(varargin{1});
ocTrj=octrajectory(trj);

if asymflag(2)
    ocTrj=ocasymptotic(ocTrj,limitset(varargin{2}));
end