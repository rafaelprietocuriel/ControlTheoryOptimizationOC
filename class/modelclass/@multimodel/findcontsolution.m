function [solObj,ocObjf,contidx]=findcontsolution(ocObj,idx,val,varargin)
%
% 

solObj=[];
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
contRes=contresult(ocObj,idx,varargin{2:end});
if iscell(contRes{1})
    contRes{1}=contRes{1}{1};
end
contSol=contRes{1}.ContinuationSolution;
contpar=zeros(1,length(contSol));

for ii=1:length(contSol)
    if isfield(contRes{1},'ExtremalSolution')
        classfield=class(contRes{1}.ExtremalSolution);
    elseif isfield(contRes{1},'LimitCycle')
        classfield=class(contRes{1}.LimitCycle);
    end
    switch classfield
        case 'octrajectory'
            hocTrj=octrajectory(contSol(ii));
            switch id
                case 'c'
                    contpar(ii)=continuationparameter(hocTrj);
                case 'i'
                    x=state(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'l'
                    x=costate(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'e'
                    x=state(ocObj,hocTrj);
                    contpar(ii)=x(val(1),end);
                case 'm'
                    par=modelparameter(hocTrj);
                    contpar(ii)=par(val(1));
            end
        case 'ocasymptotic'
            hocTrj=octrajectory(contSol(ii));
            switch id
                case 'c'
                    contpar(ii)=continuationparameter(hocTrj);
                case 'i'
                    x=state(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'l'
                    x=costate(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'e'
                    x=state(ocObj,hocTrj);
                    contpar(ii)=x(val(1),end);
                case 'm'
                    par=modelparameter(hocTrj);
                    contpar(ii)=par(val(1));
            end
        case 'mmultipath'
            
            hocTrj=sol2mmultipath(contSol(ii),ocObj);
            switch id
                case 'c'
                    contpar(ii)=continuationparameter(hocTrj(1));
                case 'i'
                    x=state(ocObj,hocTrj(1));
                    contpar(ii)=x(val(1));
                case 'm'
                    par=modelparameter(hocTrj(1));
                    contpar(ii)=par(val(1));
                case 'a'
                    arcint=arcinterval(hocTrj);
                    arcint=[arcint{:}];
                    arcint(diff(arcint)==0)=[];
                    contpar(ii)=arcint(val(1));
            end
        case 'occomposite'
            
            hocTrj=sol2occomposite(contSol(ii),ocObj);
            switch id
                case 'c'
                    contpar(ii)=continuationparameter(hocTrj(1));
                case 'i'
                    x=state(ocObj,hocTrj(1));
                    contpar(ii)=x(val(1));
                case 'm'
                    par=modelparameter(hocTrj(1));
                    contpar(ii)=par{val(3)}(val(1));
                case 'a'
                    arcint=arcinterval(hocTrj);
                    arcint=[arcint{:}];
                    arcint(diff(arcint)==0)=[];
                    contpar(ii)=arcint(val(1));
            end
        case 'dynprimitive'
            
            hocTrj=dynprimitive(contSol(ii));
            switch id
                case 'c'
                    contpar(ii)=continuationparameter(hocTrj);
                case 'i'
                    x=state(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'm'
                    par=modelparameter(hocTrj);
                    contpar(ii)=par(val(1));
            end
    end
end
contidx=cont2idx(contpar,val(2));
if isempty(contidx)
    return
end

for ii=1:length(contidx)
    if length(contidx)>1
        switch classfield
            case 'octrajectory'
                solObj{ii}=octrajectory(contSol(contidx(ii)));
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));
                
            case 'ocasymptotic'
                ocLim=contRes{1}.LimitSet;
                solObj{ii}=ocasymptotic(octrajectory(contSol(contidx(ii))),ocLim);
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));

            case 'mmultipath'
                contSol0{ii}=contSol(contidx(ii));
                solObj{ii}=sol2mmultipath(contSol(contidx(ii)),solObj);
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));

            case 'occomposite'
                contSol0{ii}=contSol(contidx(ii));
                solObj{ii}=sol2occomposite(contSol(contidx(ii)),ocObj);
                ocObjf{ii}=[];%changeparametervalue(ocObj,modelparameter(solObj{ii}));
         case 'dynprimitive'
                contSol0{ii}=contSol(contidx(ii));
                solObj{ii}=dynprimitive(contSol(contidx(ii)));
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));
            
       end
    else
        switch classfield
            case 'octrajectory'
                solObj=octrajectory(contSol(contidx(ii)));
                ocObjf=changeparametervalue(ocObj,modelparameter(solObj));
            case 'ocasymptotic'
                solObj=octrajectory(contSol(contidx(ii)));
                ocObjf=changeparametervalue(ocObj,modelparameter(solObj));
                if isfield(contRes{1},'LimitSet')
                    ocLim=contRes{1}.LimitSet;
                else
                    eq.y=contSol(contidx(ii)).solverinfo.parameters(contSol(contidx(ii)).solverinfo.equilibriumcoord);
                    eq.arcarg=arcargument(limitset(contRes{1}.ExtremalSolution));
                    ocLim=dynprimitive(eq);
                    ocLim.linearization=linearize(ocLim,ocObjf,'dependentvar',1);
                end
                solObj=ocasymptotic(solObj,ocLim);
            case 'mmultipath'

                contSol0=contSol(contidx(ii));
                solObj=sol2mmultipath(contSol(contidx(ii)),ocObj);
                ocObjf=changeparametervalue(ocObj,modelparameter(solObj));
         case 'dynprimitive'
                contSol0=contSol(contidx(ii));
                solObj=dynprimitive(contSol(contidx(ii)));
                ocObjf=changeparametervalue(ocObj,modelparameter(solObj{ii}));
            case 'occomposite'
                contSol0=contSol(contidx);
                solObj=sol2occomposite(contSol(contidx),ocObj);
                ocObjf=[];%changeparametervalue(ocObj,modelparameter(solObj{ii}));
        end
    end
end
