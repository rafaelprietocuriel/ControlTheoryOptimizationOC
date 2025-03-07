function b=implicitcontrols(ocStruct,arcid)
%
%
b=[];
if isempty(ocStruct)
    return
end
if nargin==1
    arcid=[];
end
if isempty(arcid)
    arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
    arcid=arcidentifier.value;
end
if ~iscell(arcid)
    arcid=cellstr(num2str(arcid));
end
b=false;
for ii=1:length(arcid)
    out=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcid{ii});
    if ~isempty(out.value)
        b=true;
        return
    end
end
