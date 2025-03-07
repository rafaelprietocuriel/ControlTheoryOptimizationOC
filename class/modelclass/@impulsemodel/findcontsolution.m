function [solObj,ocObjf]=findcontsolution(ocObj,idx,val,varargin)
%
% 

solObj=[];
ocObjf=[];

if isempty(ocObj)
    return
end
contRes=contresult(ocObj);
contSol=contRes{idx}.ContinuationSolution;
contpar=zeros(1,length(contSol));

for ii=1:length(contSol)
    switch class(contRes{idx}.ExtremalSolution)
        case 'hybridoctrajectory'
            hocTrj=hybridoctrajectory(contSol(ii));
            contpar(ii)=continuationparameter(hocTrj);
    end
end
contidx=cont2idx(contpar,val);
if isempty(contidx)
    return
end

for ii=1:length(contidx)
    if length(contidx)>1
        switch class(contRes{idx}.ExtremalSolution)
            case 'hybridoctrajectory'
                solObj{ii}=hybridoctrajectory(contSol(contidx(ii)));
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));
        end
    else
        switch class(contRes{idx}.ExtremalSolution)
            case 'hybridoctrajectory'
                solObj=hybridoctrajectory(contSol(contidx(ii)));
                ocObjf=changeparametervalue(ocObj,modelparameter(solObj));
        end
    end
end
