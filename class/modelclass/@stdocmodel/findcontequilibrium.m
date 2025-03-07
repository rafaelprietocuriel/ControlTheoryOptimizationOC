function [y0,ocObjf,contidx]=findcontequilibrium(ocObj,idx,val,varargin)
%
% 

y0=[];
ocObjf=[];
id='';
contidx=[];
if isempty(ocObj)
    return
end
if nargin==4
    id=varargin{1};
end
if isempty(id)
    id='c';
end
contRes=matcontresult(ocObj,idx,varargin{2:end});
contSol=contRes{1}.ContinuationSolution;
contpar=zeros(1,length(contSol));

par=contSol.modelparameter;
varyparameterindex=contSol.userinfo.varyparameterindex;
y=contSol.y;

switch contRes{1}.ContinuationClassification
    case 'modelequilibrium'
        switch id
            case 'c'
                contpar=contSol.userinfo.varyparametervalue;
        end
end
contidx=cont2idx(contpar,val);
if isempty(contidx)
    return
end

y0=y(:,contidx);
ocObjf=cell(1,length(contidx));
for ii=1:length(contidx)
    par0=par;
    par0(varyparameterindex)=contpar(contidx(ii));
    ocObjf{ii}=changeparametervalue(ocObj,par0);
end