function [vart varargout]=dependentvariable(ocObj,arcarg,varargin)
%
iterate=[];
if nargin==1
    arcarg=0;
end
if nargin>=3
    iterate=varargin{1};
end
if isnumeric(arcarg)
    arcarg=num2str(arcarg);
end
if isempty(ocObj)
    return
end
if isempty(iterate)
    iterate=1;
end
info=retrievediffmodelinformation(ocObj.Model,'equationvariablenamet',arcarg);
vart=info.value;

for ii=1:iterate
    varargout{ii}=strrep(vart,'_t',['_tp' num2str(ii)]);
end
