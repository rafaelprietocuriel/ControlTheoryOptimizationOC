function transStruct=generatevariablestruct(ocStruct,transStruct,varargin)

switch modeltype(ocStruct)
    case 'standardmodel'
         transStruct=generatevariablestruct4standardmodel(ocStruct,transStruct,varargin{:});
    case 'multistagemodel'
         transStruct=generatevariablestruct4multistagemodel(ocStruct,transStruct,varargin{:});
    case 'impulsemodel'
         transStruct=generatevariablestruct4impulsemodel(ocStruct,transStruct,varargin{:});
    case 'standarddiffmodel'
         transStruct=generatevariablestruct4standarddiffmodel(ocStruct,transStruct,varargin{:});
    case 'odemodel'
         transStruct=generatevariablestruct4odemodel(ocStruct,transStruct,varargin{:});
    case 'ppdemodel'
         transStruct=generatevariablestruct4ppdemodel(ocStruct,transStruct,varargin{:});
    case 'staticoptmodel'
         transStruct=generatevariablestruct4staticoptmodel(ocStruct,transStruct,varargin{:});
    case 'differentialgame'
         transStruct=generatevariablestruct4differentialgame(ocStruct,transStruct,varargin{:});
end