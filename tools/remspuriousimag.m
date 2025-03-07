function ocrealSolcell=remspuriousimag(ocSolcell,opt)
%
% REMSPURIOUSIMAG
if nargin==1
    opt=defaultocoptions;
end

if ~iscell(ocSolcell) 
    if ~isocsolutionclass(ocSolcell)
        ocmaterror('First argument is not a cell or OCMAT class.')
    else
        ocSolcell{1}=ocSolcell;
    end
end

ocrealSolcell=ocSolcell;
for ii=1:numel(ocSolcell) 
    if ~isocsolutionclass(ocSolcell{ii})
        ocmaterror('First argument is not a cell of solution classes.')
    end
    ocrealSolcell{ii}=remspuriousimag(ocSolcell{ii},opt);
end