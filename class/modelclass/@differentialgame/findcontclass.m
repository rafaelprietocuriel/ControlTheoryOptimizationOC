function idx=findcontclass(dgObj,contclass,varargin)
%

filename=[];
if nargin>=3
    filename=varargin{1};
end

if isempty(filename)
    contresultStruct=contresult(dgObj);
    if isempty(contresultStruct)
        idx=[];
        return
    end

    idx=zeros(1,length(contresultStruct));
    for ii=1:length(contresultStruct)
        idx(ii)=any(strcmp(contclass,contresultStruct{ii}.ContinuationClassification));
    end
    idx=find(idx);
else
    idx=cell(1,length(filename));
    for ii=1:length(filename)
        if iscell(filename)
            load(dgObj,[],1,filename{ii});
        elseif isstruct(filename) && isfield(filename,'name')
            load(dgObj,[],1,filename(ii).name);
        end
        idx{ii}=findcontclass(dgObj,contclass);
    end
end