function [idx,ocEx]=findextremal2ep(ocObj,ocEP,varargin)
%
% FINDEXTREMAL2EP returns index of extremals that converge to a specific
% equilibrium

ocEx=[];
idx=[];
contidx=[];
opt=[];
if isempty(ocObj)
    return
end
if nargin>=3
    contidx=varargin{1};
end
if nargin>=4
    opt=varargin{2};
end
contRes=contresult(ocObj);
if ~isoctrajectory(ocEP) 
    if isnumeric(ocEP) && length(ocEP)==1
        epidx=ocEP;
        ocEP=equilibrium(ocObj);
        ocEP=ocEP{epidx};
    else
        return
    end
elseif ~isequilibrium(ocEP)
    return
end
if isempty(ocEP)
    return
end
ocEx=extremalsolution(ocObj);
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(contidx)
    contidx=1:length(contRes);
end
idx=zeros(1,length(contidx));
counter=0;
for ii=contidx
    switch class(contRes{ii}.ExtremalSolution)
        case 'ocasymptotic'
            ocEPcmp{1}=ocEP;
            ocEPcmp{2}=limitset(contRes{ii}.ExtremalSolution);
            ocEPcmp=uniqueoc(ocEPcmp,opt);
            if length(ocEPcmp)==1
                counter=counter+1;
                idx(counter)=ii;
            end
        case 'ocmultipath'
    end
end
idx(idx==0)=[];
ocEx=ocEx(idx);