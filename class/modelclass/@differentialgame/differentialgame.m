function dgObj=differentialgame(firstarg,Result,varargin)
%
% DIFFERENTIALGAME initializes a standard open-loop game. This object consists of the
% two fields MODEL and RESULT. 
%
% DIFFERENTIALGAME() an empty object is returned
%
% DIFFERENTIALGAME(OCSTRUCT) 'OCSTRUCT' is a regular model structure, as e.g.
% returned by the initialization process 
%
% DIFFERENTIALGAME(MODELNAME) the model structure is loaded from the data file
% stored in the corresponding user data folder and generated during the
% initialization process.
%
% DIFFERENTIALGAME(FILENAME) the model structure is loaded from the data file
% 'FILENAME', where 'FILENAME' is either on the MATLAB path or has to
% consist of the full path.
%
% DIFFERENTIALGAME(.,RESULT) RESULT is a structure storing results for the model.

force=[];

% determine if symbolic toolbox is registered
if ~isempty(ver('symbolic'))
    if verLessThan('matlab','7.6.0')
        % maple kernel is used for symbolic calculations
        symkernel='maple';
    else
        % mupad kernel is used for symbolic calculations
        symkernel='mupad';
    end
else
    symkernel='none';
end

if nargin==0 || isempty(firstarg)
    firstarg=struct([]);
    Result=struct([]);
end
if nargin<=1
    Result=struct([]);
end
if nargin>=3
    force=varargin{1};
end
if isempty(force)
    force=0;
end
% test if ocStruct is a structure, if ocStruct is a model name load the
% corresponding model data
if ischar(firstarg)
    [ocStruct filename]=loadmodeldata(firstarg);
    if ~isempty(ocStruct)
        ocmatmsg('Model data structure successfully loaded from data file:\n%s\n',filename)
    end
elseif isstruct(firstarg)
    ocStruct=firstarg;
else
    ocmaterror('First input argument ''%s'' is neither a model structure nor a model name.',inputname(1))
end
if ~isempty(ocStruct) && ~force
    if ~isstruct(ocStruct) || differentialgametestconsistency(ocStruct,'basicmodelstructure',symkernel)
        ocmaterror('Initialization aborted. The model structure is not consistent.\n')
    end
    if isfield(ocStruct,'foc')
        ocStruct=rmfield(ocStruct,'foc');
    end
else
    %ocmatmsg('Empty model returned!\n')
end
if ~isstruct(Result) 
    %ocmatmsg('The second input argument is not a structure, set to empty.\n')
    Result=struct([]);
end
dgObj.Model=ocStruct;
dgObj.Result=Result;
dgObj=class(dgObj,'differentialgame');