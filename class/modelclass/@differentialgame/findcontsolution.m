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
                case 't'
                    arcint=arcinterval(hocTrj);
                    contpar(ii)=arcint(val(1));
                case 'i'
                    x=state(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'l'
                    x=costate(ocObj,hocTrj);
                    contpar(ii)=x(val(1));
                case 'e'
                    x=hocTrj.y;
                    %x=state(ocObj,hocTrj,1);
                    contpar(ii)=x(val(1),end);
                case 'm'
                    par=modelparameter(hocTrj);
                    if iscell(val)
                        valtmp=val;
                        parindex=parameterindex(ocObj,val{1});
                        clear val
                        val(2)=valtmp{2};
                    else
                        parindex=parameterindex(ocObj,val(1));
                    end
                    contpar(ii)=par(parindex);
                case 'T'
                    arcint=arcinterval(hocTrj(1));
                    contpar(ii)=arcint(end);
            end
        case 'ocgradtrajectory'
            ocgTrj=ocgradtrajectory(contSol(ii).extremal);
            switch id
                case 'm'
                    par=modelparameter(ocgTrj);
                    if iscell(val)
                        valtmp=val;
                        parindex=parameterindex(ocObj,val{1});
                        clear val
                        val(2)=valtmp{2};
                    else
                        parindex=parameterindex(ocObj,val(1));
                    end
                    contpar(ii)=par(parindex);
                case 'i'
                    x=state(ocgTrj);
                    contpar(ii)=x(val(1));
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
                    x=costate(ocObj,hocTrj,1);
                    contpar(ii)=x(val(1));
                case 'e'
                    x=state(ocObj,hocTrj,1);
                    contpar(ii)=x(val(1),end);
                case 'm'
                    par=modelparameter(hocTrj);
                    contpar(ii)=par(val(1));
            end
        case 'ocmultipath'
            
            hocTrj=sol2ocmultipath(contSol(ii));
            switch id
                case 'c'
                    contpar(ii)=continuationparameter(hocTrj(1));
                case 'i'
                    x=state(ocObj,hocTrj(1));
                    contpar(ii)=x(val(1));
                case 'm'
                    par=modelparameter(hocTrj(1));
                    contpar(ii)=par(val(1));
                case 'T'
                    arcint=arcinterval(hocTrj(1));
                    contpar(ii)=arcint(end);
                case 'e'
                    x=hocTrj(val(3)).y;
                    contpar(ii)=x(val(1),end);
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
            case 'ocgradtrajectory'
                solObj{ii}=ocgradtrajectory(contSol(contidx(ii)).extremal);
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));
                
            case 'ocasymptotic'
                solObj{ii}=octrajectory(contSol(contidx(ii)));
                ocObjf{ii}=changeparametervalue(ocObj,modelparameter(solObj{ii}));
                if isfield(contRes{1},'LimitSet')
                    ocLim=contRes{1}.LimitSet;
                else
                    eq.y=contSol(contidx(ii)).solverinfo.parameters(contSol(contidx(ii)).solverinfo.equilibriumcoord);
                    eq.arcarg=arcargument(limitset(contRes{1}.ExtremalSolution));
                    ocLim=dynprimitive(eq);
                    ocLim.linearization=linearize(ocLim,ocObjf{ii},'dependentvar',1);
                end
                solObj{ii}=ocasymptotic(octrajectory(contSol(contidx(ii))),ocLim);

            case 'ocmultipath'
                for jj=1:multiplicity(contRes{1}.ExtremalSolution)
                    switch  class(contRes{1}.ExtremalSolution(jj))
                        case 'ocasymptotic'
                            ocLim{jj}=limitset(contRes{1}.ExtremalSolution(jj));
                        case 'octrajectory'
                            ocLim=[];
                    end
                end
                contSol0{ii}=contSol(contidx(ii));
                solObj{ii}=sol2ocmultipath(contSol(contidx(ii)),[],ocLim);
                par=modelparameter(solObj{ii});
                ocObjf{ii}=changeparametervalue(ocObj,par{1});
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
            case 'ocgradtrajectory'
                solObj=ocgradtrajectory(contSol(contidx).extremal);
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
            case 'ocmultipath'

                contSol0=contSol(contidx(ii));
                solObj=sol2ocmultipath(contSol(contidx(ii)),[],[]);
                par=modelparameter(solObj);
                ocObjf=changeparametervalue(ocObj,par{1});
                for jj=1:multiplicity(contRes{1}.ExtremalSolution)
                    switch  class(contRes{1}.ExtremalSolution(jj))
                        case 'ocasymptotic'
                            if isfield(contSol(contidx(ii)).solverinfo,'equilibriumcoord')
                                eq.y=contSol(contidx(ii)).solverinfo.parameters(contSol(contidx(ii)).solverinfo.equilibriumcoord{jj});
                                eq.arcarg=arcargument(limitset(contRes{1}.ExtremalSolution(jj)));
                                ocLim{jj}=dynprimitive(eq);
                                ocLim{jj}.linearization=linearize(ocLim{jj},ocObjf,'dependentvar',1);
                            else
                                ocLim{jj}=limitset(contRes{1}.ExtremalSolution(jj));
                            end
                        case 'octrajectory'
                            ocLim=[];
                    end
                end
                solObj=sol2ocmultipath(contSol(contidx(ii)),[],ocLim);
         case 'dynprimitive'
                contSol0=contSol(contidx);
                solObj=dynprimitive(contSol(contidx));
                ocObjf=changeparametervalue(ocObj,modelparameter(solObj));
        end
    end
end
