function out=arcnum(ocComp,uniformoutputflag)

if nargin==1
    uniformoutputflag=0;
end
out=[];
if isempty(ocComp)
    return
end
out=cell(1,ocComp.order);
for ii=1:ocComp.order
    out{ii}=arcnum(ocComp.path{ii});
end

if uniformoutputflag
    out=[out{:}];
end