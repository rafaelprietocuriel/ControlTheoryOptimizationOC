function idx=findcontclass(ocObj,contclass,varargin)
%

filename=[];
if nargin>=3
    filename=varargin{1};
end

if isempty(filename)
contresultStruct=contresult(ocObj);
if isempty(contresultStruct)
    idx=[];
    return
end

idx=zeros(1,length(contresultStruct));
for ii=1:length(contresultStruct)
    idx(ii)=strcmp(contclass,contresultStruct{ii}.ContinuationClassification);
end
idx=find(idx);
else
    idx=cell(1,length(filename));
    for ii=1:length(filename)
        if iscell(filename)
            load(ocObj,[],1,filename{ii});
        elseif isstruct(filename) && isfield(filename,'name')
            load(ocObj,[],1,filename(ii).name);
        end
        idx{ii}=findcontclass(ocObj,contclass);
    end
end