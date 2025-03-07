function transStruct=makestandardvariable(ocStruct,varStruct,transStruct,varargin)

switch modeltype(ocStruct)
    case 'standardmodel'
         transStruct=makestandardvariable4standardmodel(varStruct,transStruct,varargin{:});
    case 'multistagemodel'
         transStruct=makestandardvariable4multistagemodel(varStruct,transStruct,varargin{:});
    case 'differentialgame'
         transStruct=makestandardvariable4differentialgame(varStruct,transStruct,varargin{:});
    case 'odemodel'
         transStruct=makestandardvariable4standardmodel(varStruct,transStruct,varargin{:});
end