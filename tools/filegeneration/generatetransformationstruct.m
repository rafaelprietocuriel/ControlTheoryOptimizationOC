function transStruct=generatetransformationstruct(ocStruct,transStruct,varargin)

switch modeltype(ocStruct)
    case 'standardmodel'
        transStruct=generatetransformationstruct4standardmodel(ocStruct,transStruct,varargin{:});
    case 'multistagemodel'
        transStruct=generatetransformationstruct4multistagemodel(ocStruct,transStruct,varargin{:});
    case 'impulsemodel'
        transStruct=generatetransformationstruct4impulsemodel(ocStruct,transStruct,varargin{:});
    case 'standarddiffmodel'
        transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,varargin{:});
    case 'odemodel'
        transStruct=generatetransformationstruct4odemodel(ocStruct,transStruct,varargin{:});
    case 'ppdemodel'
        transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,varargin{:});
    case 'staticoptmodel'
        transStruct=generatetransformationstruct4staticoptmodel(ocStruct,transStruct,varargin{:});
    case 'differentialgame'
        transStruct=generatetransformationstruct4differentialgame(ocStruct,transStruct,varargin{:});
end