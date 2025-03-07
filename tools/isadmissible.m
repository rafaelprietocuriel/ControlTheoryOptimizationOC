function [b,adminfo]=isadmissible(ocSol,ocObj,varargin)
%
% ISADMISSIBLE calls the isadmissible function for the specific solution
% classes
% the input argument 'ocSolcell' has to be a cell of oc solution classes,
% e.g., dynprimitve, ocasymptotic, ...
b=[];
%adminfo=struct([]);
if isempty(ocObj)
    return
end

if ~iscell(ocSol) 
    if ~isocsolutionclass(ocSol)
        ocmaterror('First argument is not a cell or OCMAT class.')
    else
        ocSolcell{1}=ocSol;
    end
else
    ocSolcell=ocSol;
end

b=zeros(1,length(ocSolcell));
for ii=1:numel(ocSolcell) 
    if ~isocsolutionclass(ocSolcell{ii})
        ocmaterror('First argument is not a cell of solution classes.')
    end
    [b(ii),adminfo(ii)]=isadmissible(ocSolcell{ii},ocObj,varargin{:});
end