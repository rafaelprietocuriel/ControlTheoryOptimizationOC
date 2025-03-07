function mmObj=multimodel(firstarg,modelname,varargin)
%
% MULTIMODEL initializes a class combining multiple models. This object 
% consists of the two fields MODEL and RESULT. 
%
% MULTIMODEL() an empty object is returned
%
% MULTIMODEL(OCSTRUCT) 'OCSTRUCT' is a cell of model structures, as e.g.
% returned by the initialization process 
%
% MULTIMODEL(MODELNAME) the model structure is loaded from the data file
% stored in the corresponding user data folder and generated during the
% initialization process.
%
% MULTIMODEL(FILENAME) the model structure is loaded from the data file
% 'FILENAME', where 'FILENAME' is either on the MATLAB path or has to
% consist of the full path.
%
% MULTIMODEL(.,RESULT) RESULT is a structure storing results for the model.


if nargin==0 || isempty(firstarg)
    firstarg=struct([]);
end
% test if ocStruct is a structure, if ocStruct is a model name load the
% corresponding model data
if isstruct(firstarg)
    firstarg=struct2cell(firstarg);
end
if ~iscell(firstarg)
    ocmatmsg('First argument is not a cell or not empty\n')
    mmObj=multimodel([]);
    return
end
ocStruct=cell(1,length(firstarg));
filename=ocStruct;
NumberOfModels=length(firstarg);
for ii=1:NumberOfModels
    if ischar(firstarg{ii})
        [ocStruct{ii},filename{ii}]=loadmodeldata(firstarg{ii});
        if ~isempty(ocStruct{ii})
            ocmatmsg('Model data structure successfully loaded from data file:\n%s\n',filename{ii})
        end
    elseif isstruct(firstarg{ii})
        ocStruct{ii}=firstarg{ii};
    else
        ocmaterror('First input argument ''%s'' is neither a model structure nor a model name.',inputname(1))
    end
end
for ii=1:NumberOfModels
    if ~isempty(ocStruct{ii})
        if ~isstruct(ocStruct{ii})% || standardmodeltestconsistency(ocStruct,'basicmodelstructure',symkernel)
            ocmaterror('Initialization aborted. The model structure is not consistent.\n')
        end
        if isfield(ocStruct{ii},'foc')
            ocStruct{ii}=rmfield(ocStruct{ii},'foc');
        end
    else
        %ocmatmsg('Empty model returned!\n')
    end
end
for ii=1:NumberOfModels
    switch ocStruct{ii}.modeltype
        case 'standardmodel'
            mmObj.Model{ii}=stdocmodel(ocStruct{ii});
        case 'multistagemodel'
            mmObj.Model{ii}=multistagemodel(ocStruct{ii},ii);
        otherwise

    end
end
mmObj.NumberOfModels=NumberOfModels;
mmObj.Result=struct([]);
mmObj.modelname=modelname;
mmObj=class(mmObj,'multimodel');