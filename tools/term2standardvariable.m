function transterm=term2standardvariable(gocObj,term,varargin)%varStruct,transStruct,transformtype,varargin)
%
% transformtype: 0 ... symbolic, 1 ... array, 2 ... vectorized
% multline: if string is a matrix already decomposed in multiline format
% string is a cell array of strings
% if string contains ';' it is broken up into cells (standard matrix
% format)
% arcidentifier: is a string containing a digit identifying the specific
% arc. If arcidientifer is not empty the transformed variable names for
% this arc used for replacing. If arcidentifier is empty the standard
% transformation names are used for replacing

arcidentifier=[];
transformtype=[];
arcdependent=[];
vectorizeflag=[];
if nargin>=3
    arcidentifier=varargin{1};
end
if nargin>=4
    transformtype=varargin{2};
end
if nargin>=5
    vectorizeflag=varargin{3};
end
if nargin>=6
    arcdependent=varargin{4};
end

if isempty(arcidentifier)
    arcidentifier=0;
end
if isempty(transformtype)
    transformtype=1;
end
if isempty(vectorizeflag)
    vectorizeflag=0;
end
if isempty(arcdependent)
    arcdependent=0;
end
if isstruct(gocObj)
    ocStruct=gocObj;
    modeltype='stdocmodel';
else
    load(gocObj,'modeldata')
    modeltype=class(gocObj);
end
switch modeltype
    case 'stdocmodel'
        transStruct=generatetransformationstruct4standardmodel(ocStruct,[],'composedvar');
    case 'impulsemodel'
        transStruct=generatetransformationstruct4impulsemodel(ocStruct,[],'composedvar');
    case 'odemodel'
        transStruct=generatetransformationstruct4odemodel(ocStruct,[],'composedvar');
    otherwise
        return
end
%
varStruct.term.type='algterm';
varStruct.term.vectorize=0;
if ~arcdependent
    varStruct.term.arcidentifier=num2str(arcidentifier);
end
varStruct.term.arcdependent=arcdependent;
%
switch class(term)
    case 'char'
        if size(term,1)==1
            string{1}=term;
        else
            for ii=1:length(term)
                string{ii}=term(ii,:);
            end
        end
    case 'sym'
        for ii=1:size(term,1)
            string{ii}=removematrixstring(char(term(ii,:)));
        end
    case 'cell'
        string=term;
end
if vectorizeflag
    for ii=1:length(string)
        string{ii}=vectorize(string{ii});
    end
end
if arcdependent
    for ii=1:length(string)
        varStruct.term.arcidentifier{ii}{1}=num2str(arcidentifier);
    end
end

varStruct.term.string=string;

tvarStruct=makestandardvariable(ocStruct,varStruct,transStruct,transformtype);
transterm=tvarStruct.term.string;