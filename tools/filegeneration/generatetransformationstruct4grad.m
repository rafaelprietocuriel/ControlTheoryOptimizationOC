function transStruct=generatetransformationstruct4grad(ocStruct,transStruct,varargin)

switch modeltype(ocStruct)
    case 'standardmodel'
        transStruct=generatetransformationstruct4standardmodelgrad(ocStruct,transStruct,varargin{:});
end