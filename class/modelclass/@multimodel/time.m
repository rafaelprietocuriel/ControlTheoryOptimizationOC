function out=time(mmObj,solObj,varargin)

out=[];
connectflag=[];
if isempty(mmObj)
    return
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end

if ~ismmultipath(solObj)
    return
end
if connectflag
    for ii=1:numberofparts(solObj)
        out=[out time(mmObj.Model{ii},solObj(ii),connectflag)];
    end
else
    out=cell(1,numberofparts(solObj));
    for ii=1:numberofparts(solObj)
        out{ii}=time(mmObj.Model{ii},solObj(ii),1);
    end
end