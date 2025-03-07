function varargout=store(ocObj,varargin)
%
% STORE results of continuation process
%
% STORE(OCOBJ,PROBLEMTYPE) store results of continuation process
% given by PROBLEMTYPE:
%   'extremal2ep'   ... saddle-path of an equilibrium, continuing along the
%                       initial point
%   'extremalp2ep'  ... saddle-path of an equilibrium, continuing along a
%                       parameter value
%   'extremalt2ep'  ... saddle-path of an equilibrium, continuing the
%                       truncation time
%   'indifferencesolution'  ... continuation of an indifference threshold
%                       (Skiba point)
%   'limitextremal' ... continuation of a limitpoint solution
%
% STORE(OCOBJ,OCELEMENT) stores the object OCELEMENT in the Result field of
% the stdocmodel object OCOBJ.
% Possible classes for OCELEMENT are
%        dynprimitive: equilibrium or limit cycle
%        ocasymptotic: asymptotic solutions
%        cell array of 'dynprimitives' or 'ocasymptotics'
%
% The default names of the Result field are
%       Result.Equilibrium: for a dynprimitive being an equilibrium
%       Result.LimitCycle: for a dynprimitive being a limit cycle
%       Result.ExtremalSolution: for an ocasymptotic
%       Result.Continuation: for the result of a continuation process
% If no output argument is specified it is written into the object OCOBJ.
%
% OCOBJN=STORE(...) the result is returned to the new stdocmodel instance
% OCOBJN


if isempty(ocObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1}) || isoctrajectory(varargin{1}) || isoccurve(varargin{1}) ||  isocmultipath(varargin{1}) || isocgradlimitset(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(ocObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),ocObj);
        else
            varargout{1}=ocObj;
        end
        return
    else
        ocmaterror('Input argument ''%s'' needs the specification of a field name.\n',inputname(2))
    end
    fieldvalue=varargin{1};
end

if ischar(varargin{1}) && nargin==2
    fieldname=varargin{1};
    fieldvalue=[];
end

if ischar(varargin{1}) && nargin==3
    fieldname=varargin{1};
    fieldvalue=varargin{2};
end

if isempty(fieldname) && isempty(fieldvalue)
    ocmaterror('Wrong input arguments')
end
[ocStruct,ocResultFieldName]=generateelement(ocObj,fieldname,fieldvalue);
if isempty(ocResultFieldName) || isempty(ocStruct)
    return
end
if isfield(ocObj.Result,ocResultFieldName)
    ocObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    ocObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),ocObj);
else
    varargout{1}=ocObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(ocObj,fieldname,fieldvalue)
ocResultStruct=[];
ocResultFieldName='';
switch fieldname
    case {'dynprimitive','gdynprimitive'}
        if isequilibrium(fieldvalue)
            ocResultStruct=fieldvalue;
            ocResultFieldName='Equilibrium';
        else
            ocResultStruct=fieldvalue;
            ocResultFieldName='PeriodicSolution';
        end
    case 'ocgradlimitset'
        if isequilibrium(fieldvalue)
            ocResultStruct=fieldvalue;
            ocResultFieldName='GradEquilibrium';
        else
            ocResultStruct=fieldvalue;
            ocResultFieldName='GradPeriodicSolution';
        end
    case 'ocasymptotic'
        ocResultStruct=fieldvalue;
        ocResultFieldName='ExtremalSolution';

    case {'extremal2ep','extremalt2ep','extremale2ep','extremalc2ep','extremal2emf','extremalt2emf','extremal2inf','extremalt2inf','saddlepath2ep'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        mn=modelname(ocObj);

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        solver=sout(1).data.sol.solver;
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(1:OCBVP.numode(1),counter)=ocExStruct(ii).y(1:OCBVP.numode(1),1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                if OCMATAE.objectivevaluecalc
                    slMfStruct.userinfo.objectivevalue(counter)=ocExStruct(ii).y(OCMATAE.objectivevaluecoord,end);
                end
            end
        end
        
        slMfStruct.arcposition=[];
        if strcmp(fieldname,'extremale2ep')
            ocEP=dynprimitive(OCMATAE.saddlepoint,OCMATAE.endarcarg,OCMATAE.linearization);
        elseif any(strcmp(fieldname,{'extremalt2ep','extremalc2ep','saddlepath2ep'}))
            ocEP=dynprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.linearization);
        elseif any(strcmp(fieldname,{'extremal2inf','extremalt2inf'}))
            if ~isempty(OCMATAE.EP.saddlepoint)
                ocEP=generateequilibrium(OCMATAE.EP.saddlepoint,OCMATCONT.HE.arcarg(end),parametervalue(ocObj),mn,'dynprimitive');
            else
                ocEP=dynprimitive();
            end
        elseif strcmp(fieldname,'extremal2ep')
            switch OCMATCONT.continuationtype
                case  {'initialstate','time'}
                    switch solver
                        case 'gbvp4c'
                            ocEP=generateequilibrium(OCMATAE.EP.saddlepoint,OCMATAE.EP.arcarg,parametervalue(ocObj),mn,'gdynprimitive');
                            ocEPN=ocEP;
                        otherwise
                            ocEP=generateequilibrium(OCMATAE.EP.saddlepoint,OCMATAE.EP.arcarg,parametervalue(ocObj),mn,'dynprimitive');
                            ocEPN=ocEP;
                    end
                otherwise
                    ocEP=[];
            end
        end
        if ~isempty(soln)
            ocObjTmp=ocObj;
            if any(strcmp(fieldname,{'extremal2emf','extremalt2emf'}))
                ocEPN=generateequilibrium(soln.solverinfo.saddlepoint,OCMATCONT.HE.arcarg(end),soln.modelparameter,mn,'dynprimitive');
                ocEP=generateequilibrium(sol.solverinfo.saddlepoint,OCMATCONT.HE.arcarg(end),soln.modelparameter,mn,'dynprimitive');
            elseif strcmp(fieldname,'extremal2inf')
                if strcmp(OCMATCONT.continuationtype,'parameter')
                else
                    ocEPN=ocEP;
                end
            elseif strcmp(fieldname,'extremal2ep')
                if strcmp(OCMATCONT.continuationtype,'parameter')
                    switch solver
                        case 'gbvp4c'
                            ocEPN=generateequilibrium(soln.solverinfo.saddlepoint,OCMATAE.EP.arcarg,soln.modelparameter,mn,'gdynprimitive');
                            ocEP=generateequilibrium(sol.solverinfo.saddlepoint,OCMATAE.EP.arcarg,sol.modelparameter,mn,'gdynprimitive');
                        otherwise
                            ocEPN=generateequilibrium(soln.solverinfo.saddlepoint,OCMATAE.EP.arcarg,soln.modelparameter,mn,'dynprimitive');
                            ocEP=generateequilibrium(sol.solverinfo.saddlepoint,OCMATAE.EP.arcarg,sol.modelparameter,mn,'dynprimitive');
                    end
                end
            end
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgasymptotic(ocgtrajectory(soln,soln.solverinfo.odenum,ocObjTmp),ocEPN);
                otherwise
                    ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,ocObjTmp),ocEPN);
            end

            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgasymptotic(ocgtrajectory(sol,sol.solverinfo.odenum,ocObjTmp),ocEP);
                otherwise
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObjTmp),ocEP);
            end
            
        else
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        if ~isempty(soltarget)
            ocObjTmp=ocObj;
            if any(strcmp(fieldname,{'extremal2emf','extremalt2emf'}))
                ocEP=generateequilibrium(soltarget.solverinfo.saddlepoint,OCMATCONT.HE.arcarg(end),soltarget.modelparameter,soltarget.modelname,'dynprimitive');
            elseif strcmp(fieldname,'extremal2ep')
                if strcmp(OCMATCONT.continuationtype,'parameter')
                    switch solver
                        case 'gbvp4c'
                            ocEP=generateequilibrium(soltarget.solverinfo.saddlepoint,OCMATAE.EP.arcarg,soltarget.modelparameter,soltarget.modelname,'gdynprimitive');
                        otherwise
                            ocEP=generateequilibrium(soltarget.solverinfo.saddlepoint,OCMATAE.EP.arcarg,soltarget.modelparameter,soltarget.modelname,'dynprimitive');
                    end
                end
            end
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgasymptotic(ocgtrajectory(soltarget,soltarget.solverinfo.odenum,ocObjTmp),ocEP);
                otherwise
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,ocObjTmp),ocEP);
            end
            
        else
            ocObjTmp=ocObj;
            if any(strcmp(fieldname,{'extremal2emf','extremalt2emf'}))
                 ocEP=generateequilibrium(ocExStruct(end).solverinfo.saddlepoint,OCMATCONT.HE.arcarg(end),ocExStruct(end).modelparameter,ocExStruct(end).modelname,'dynprimitive');
            elseif strcmp(fieldname,'extremal2ep')
                if strcmp(OCMATCONT.continuationtype,'parameter')
                    switch solver
                        case 'gbvp4c'
                            ocEP=generateequilibrium(ocExStruct(end).solverinfo.saddlepoint,OCMATAE.EP.arcarg,ocExStruct(end).modelparameter,ocExStruct(end).modelname,'gdynprimitive');
                        otherwise
                            ocEP=generateequilibrium(ocExStruct(end).solverinfo.saddlepoint,OCMATAE.EP.arcarg,ocExStruct(end).modelparameter,ocExStruct(end).modelname,'dynprimitive');
                    end
                end
            end
            %DataAdaptation(ocExStruct(end).solverinfo.tmesh);
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgasymptotic(ocgtrajectory(ocExStruct(end),ocExStruct(end).solverinfo.odenum,ocObjTmp),ocEP);
                otherwise
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),ocObjTmp),ocEP);
            end
        end
        if any(strcmp(fieldname,{'extremal2ep','extremal2inf','extremalt2inf','extremale2ep','extremalc2ep','extremal2emf','extremalt2emf','saddlepath2ep'}))
            ocResultStruct.SliceManifold=occurve(slMfStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'extremalp2ep','extremalp2epuser','extremalp2inf'}
        global OCMATCONT OCMATAE OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        meshadaptationflag=OCMATCONT.meshadaptationflag;
        OCMATCONT.meshadaptationflag=1;
        for ii=1:numel(bvpout)
            DataAdaptation(bvpout(ii).tmesh);
            try
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        OCMATCONT.meshadaptationflag=meshadaptationflag;

        if ~isempty(soln)
            ocObjtmp=changeparametervalue(ocObj,OCMATAE.varyparameterindex,soln.solverinfo.continuationparameter);
            if any(strcmp(fieldname,{'extremalp2inf'}))
                ocEP=dynprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.linearization);
            else
                equcoord=soln.solverinfo.parametercoord(soln.solverinfo.equilibriumcoord);
                hatx.y=soln.solverinfo.coeff(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                hatx.linearization=jacobian(ocObjtmp,hatx);
                ocEP=dynprimitive(hatx);
                ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,ocObjtmp),ocEP);
                equcoord=sol.solverinfo.parametercoord(sol.solverinfo.equilibriumcoord);
                ocObjtmp=changeparametervalue(ocObj,OCMATAE.varyparameterindex,sol.solverinfo.continuationparameter);
                hatx.y=sol.solverinfo.coeff(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                if ~isempty(OCMATAE.fixcoordinate)
                    hatx.y(OCMATAE.varcoordinate)=hatx.y;
                    hatx.y(OCMATAE.fixcoordinate)=OCMATAE.fixcoordinatevalue;
                end
                hatx.linearization=jacobian(ocObjtmp,hatx);
                ocEP=dynprimitive(hatx);
            end
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObjtmp),ocEP);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                sol=soltarget;
            else
                sol=ocExStruct(end);
            end
            ocObjtmp=changeparametervalue(ocObj,OCMATAE.varyparameterindex,sol.solverinfo.continuationparameter);
            if any(strcmp(fieldname,{'extremalp2inf'}))
                ocEP=dynprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.linearization);
            else
                equcoord=sol.solverinfo.parametercoord(sol.solverinfo.equilibriumcoord);
                hatx.y=sol.solverinfo.coeff(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                if ~isempty(OCMATAE.fixcoordinate)
                    hatx.y(OCMATAE.varcoordinate)=hatx.y;
                    hatx.y(OCMATAE.fixcoordinate)=OCMATAE.fixcoordinatevalue;
                end
                hatx.linearization=jacobian(ocObjtmp,hatx);
                ocEP=dynprimitive(hatx);
            end
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObjtmp),ocEP);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
            matlabsol=transform2nativematlab(ocExStruct(end).solverinfo.tmesh,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        end
        ocResultStruct.ContinuationParameter=OCMATAE.varyparameterindex;
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);



    case {'indifferencesolution4lcae'}
        global OCMATCONT OCMATINDIF OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        meshadaptationflag=OCMATCONT.meshadaptationflag;
        OCMATCONT.meshadaptationflag=1;
        for ii=1:numel(bvpout)
            DataAdaptation(bvpout(ii).tmesh);
            try
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        OCMATCONT.meshadaptationflag=meshadaptationflag;

        if ~isempty(soln)
            sol=soln;
            arcposition=find(diff(sol.x)==0);
            leftarcindex=[1 arcposition+1];
            rightarcindex=[arcposition numel(sol.x)];

            ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate);
            hatx.x=0;
            hatx.arcarg=OCMATINDIF.arcarg{1}(end);
            hatx.linearization=jacobian(ocObjtmp,hatx);
            ocEP=dynprimitive(hatx);

            trjidx=find(OCMATINDIF.solutionindex==1);
            trjleftarcindex=leftarcindex(trjidx);
            trjrightarcindex=rightarcindex(trjidx);
            lcidx=find(OCMATINDIF.solutionindex==2);
            lcleftarcindex=leftarcindex(lcidx);
            lcrightarcindex=rightarcindex(lcidx);

            ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
            ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
            ocTrj.arcarg=OCMATINDIF.arcarg{1};
            ocTrj.arcinterval=sol.solverinfo.arcinterval{1};
            ocTrj.timehorizon=sol.solverinfo.timehorizon(1);
            ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
            ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
            ocTrj.modelparameter=sol.modelparameter;
            ocTrj.modelname=sol.modelname;

            ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
            ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
            ocLC.octrajectory.arcarg=OCMATINDIF.arcarg{2};
            ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{2};
            ocLC.octrajectory.timehorizon=sol.solverinfo.timehorizon(2);
            ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
            ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
            ocLC.octrajectory.modelparameter=sol.modelparameter;
            ocLC.octrajectory.modelname=sol.modelname;
            ocLC.period=sol.solverinfo.timehorizon(2);

            ocResultStruct.NonadmissibleLimitCycle=dynprimitive(ocLC);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);

            sol=sollast;
            arcposition=find(diff(sol.x)==0);
            leftarcindex=[1 arcposition+1];
            rightarcindex=[arcposition numel(sol.x)];

            ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate);
            hatx.x=0;
            hatx.arcarg=OCMATINDIF.arcarg{1}(end);
            hatx.linearization=jacobian(ocObjtmp,hatx);
            ocEP=dynprimitive(hatx);

            trjidx=find(OCMATINDIF.solutionindex==1);
            trjleftarcindex=leftarcindex(trjidx);
            trjrightarcindex=rightarcindex(trjidx);
            lcidx=find(OCMATINDIF.solutionindex==2);
            lcleftarcindex=leftarcindex(lcidx);
            lcrightarcindex=rightarcindex(lcidx);

            ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end))-sol.x(trjleftarcindex(1));
            ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
            ocTrj.arcarg=OCMATINDIF.arcarg{1};
            ocTrj.arcinterval=sol.solverinfo.arcinterval{1};
            ocTrj.timehorizon=sol.solverinfo.timehorizon(1);
            ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
            ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
            ocTrj.modelparameter=sol.modelparameter;
            ocTrj.modelname=sol.modelname;
            
            ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
            ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
            ocLC.octrajectory.arcarg=OCMATINDIF.arcarg{2};
            ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{2};
            ocLC.octrajectory.timehorizon=sol.solverinfo.timehorizon(2);
            ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
            ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
            ocLC.octrajectory.modelparameter=sol.modelparameter;
            ocLC.octrajectory.modelname=sol.modelname;
            ocLC.period=sol.solverinfo.timehorizon(2);
            
            ocResultStruct.LimitCycle=dynprimitive(ocLC);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocTrj,ocObjtmp),ocEP);
        
        else
            if ~isempty(soltarget)
                sol=soltarget;
            else
                sol=ocExStruct(end);
            end
            arcposition=find(diff(sol.x)==0);
            leftarcindex=[1 arcposition+1];
            rightarcindex=[arcposition numel(sol.x)];

            ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoordinate);
            hatx.x=0;
            hatx.arcarg=OCMATINDIF.arcarg{1}(end);
            hatx.linearization=jacobian(ocObjtmp,hatx);
            ocEP=dynprimitive(hatx);

            trjidx=find(OCMATINDIF.solutionindex==1);
            trjleftarcindex=leftarcindex(trjidx);
            trjrightarcindex=rightarcindex(trjidx);
            lcidx=find(OCMATINDIF.solutionindex==2);
            lcleftarcindex=leftarcindex(lcidx);
            lcrightarcindex=rightarcindex(lcidx);

            ocTrj.x=sol.x(trjleftarcindex(1):trjrightarcindex(end));
            ocTrj.y=sol.y(:,trjleftarcindex(1):trjrightarcindex(end));
            ocTrj.arcarg=OCMATINDIF.arcarg{1};
            ocTrj.arcinterval=sol.solverinfo.arcinterval{1};
            ocTrj.timehorizon=sol.solverinfo.timehorizon(1);
            ocTrj.arcposition=[trjleftarcindex;trjrightarcindex];
            ocTrj.arcposition=ocTrj.arcposition-ocTrj.arcposition(1,1)+1;
            ocTrj.modelparameter=sol.modelparameter;
            ocTrj.modelname=sol.modelname;

            ocLC.octrajectory.x=sol.x(lcleftarcindex(1):lcrightarcindex(end))-sol.x(lcleftarcindex(1));
            ocLC.octrajectory.y=sol.y(:,lcleftarcindex(1):lcrightarcindex(end));
            ocLC.octrajectory.arcarg=OCMATINDIF.arcarg{2};
            ocLC.octrajectory.arcinterval=sol.solverinfo.arcinterval{2};
            ocLC.octrajectory.timehorizon=sol.solverinfo.timehorizon(2);
            ocLC.octrajectory.arcposition=[lcleftarcindex;lcrightarcindex];
            ocLC.octrajectory.arcposition=ocLC.octrajectory.arcposition-ocLC.octrajectory.arcposition(1,1)+1;
            ocLC.octrajectory.modelparameter=sol.modelparameter;
            ocLC.octrajectory.modelname=sol.modelname;
            ocLC.period=sol.solverinfo.timehorizon(2);
            
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocTrj,ocObjtmp),ocEP);
            ocResultStruct.LimitCycle=dynprimitive(ocLC);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);

    case {'heteroclinicep2lc'}
        global OCMATCONT OCMATHET OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATHETORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        meshadaptationflag=OCMATCONT.meshadaptationflag;
        OCMATCONT.meshadaptationflag=1;
        for ii=1:numel(bvpout)
            DataAdaptation(bvpout(ii).tmesh);
            try
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        OCMATCONT.meshadaptationflag=meshadaptationflag;

        if ~isempty(soltarget)
            ocResultStruct.ExtremalSolution=heteroclinicep2lc2ocmultipath(ocObj,soltarget);
            ocResultStruct.NonadmissibleSolution=ocmultipath();
        elseif ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=heteroclinicep2lc2ocmultipath(ocObj,soln);
            ocResultStruct.ExtremalSolution=heteroclinicep2lc2ocmultipath(ocObj,sol);
        elseif ~isempty(sollast)
            ocResultStruct.ExtremalSolution=heteroclinicep2lc2ocmultipath(ocObj,sollast);
            ocResultStruct.NonadmissibleSolution=ocmultipath();
        end
        
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATHETORIGINAL,OCBVPORIGINAL);


    case {'extremal4ft','extremalt4ft','extremalp4ft','extremaltF4ft','extremalpF4ft','extremalF4ft','extremalt4inft'}
        global OCMATCONT OCMATFTE OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATFTEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        solver=sout(1).data.sol.solver;
        counter=0;
        %meshadaptationflag=OCMATCONT.meshadaptationflag;
        OCMATCONT.meshadaptationflag=1;
        maxnumode=max(OCBVP.numode);
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                ocExStruct(counter).solverinfo.stepwidth=bvpout(ii).stepwidth;
                slMfStruct.y(maxnumode,counter)=ocExStruct(ii).y(maxnumode,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                endpointStruct.y(maxnumode,counter)=ocExStruct(ii).y(maxnumode,end);
                endpointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(end);
                endpointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                if OCMATFTE.objectivevaluecalc
                    slMfStruct.userinfo.objectivevalue(counter)=ocExStruct(ii).y(OCMATFTE.objectivevaluecoord,end);
                end
            end
        end
        %OCMATCONT.meshadaptationflag=meshadaptationflag;
        slMfStruct.arcposition=[];
        endpointStruct.arcposition=[];
        if ~isempty(soltarget)
            sol=soltarget;
            if strcmp(fieldname,'extremalp4ft')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgtrajectory(sol,sol.solverinfo.odenum,ocObjtmp);
                    ocResultStruct.NonadmissibleSolution=ocgtrajectory([]);
                otherwise
                    ocResultStruct.ExtremalSolution=octrajectory(sol,ocObjtmp);
                    ocResultStruct.NonadmissibleSolution=octrajectory([]);
            end
        elseif ~isempty(soln)
            if strcmp(fieldname,'extremalp4ft')
                ocObjtmp=changeparametervalue(ocObj,soln.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgtrajectory(soln,soln.solverinfo.odenum,ocObjtmp);
                otherwise
                    ocResultStruct.NonadmissibleSolution=octrajectory(soln,ocObjtmp);
            end
            
            if strcmp(fieldname,'extremalp4ft')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgtrajectory(sol,sol.solverinfo.odenum,ocObjtmp);
                otherwise
                    ocResultStruct.ExtremalSolution=octrajectory(sol,ocObjtmp);
            end
        else
            sol=ocExStruct(end);
            if strcmp(fieldname,'extremalp4ft')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            switch solver
                case 'gbvp4c'
                    ocResultStruct.ExtremalSolution=ocgtrajectory(sol,sol.solverinfo.odenum,ocObjtmp);
                    ocResultStruct.NonadmissibleSolution=ocgtrajectory([]);
                otherwise
                    ocResultStruct.ExtremalSolution=octrajectory(sol,ocObjtmp);
                    ocResultStruct.NonadmissibleSolution=octrajectory([]);
            end
        end
        ocResultStruct.SliceManifold=occurve(slMfStruct);
        ocResultStruct.EndpointManifold=occurve(endpointStruct);
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if strcmp(fieldname,'extremalp4ft')
            ocResultStruct.ContinuationParameter=OCMATFTE.continuationindex;
            ocResultStruct.LastModel=ocObjtmp;
        end
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATFTEORIGINAL,OCBVPORIGINAL);


    case {'optisocline'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.isoclinecoordinate=OCMATAE.isoclinecoordinate;
            end
        end
        slMfStruct.arcposition=[];
        ocEP=dynprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.linearization);
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,ocObj),ocEP);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObj),ocEP);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        if ~isempty(soltarget)
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,ocObj),ocEP);
            matlabsol=transform2nativematlab(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            DataAdaptation(ocExStruct(end).solverinfo.tmesh);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),ocObj),ocEP);
            matlabsol=transform2nativematlab(ocExStruct(end).solverinfo.tmesh,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        end
        ocResultStruct.Isocline=occurve(slMfStruct);
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
    case {'modelequilibrium','modelequilibriummf','modelequilibriumu2c','modelequilibriumc'}
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            modelfolder=getocmatfolder('userdata',modeltype(ocObj),modelname(ocObj));
            OCMATFINITCONT.basicglobalvarfilename=fullocmatfile(modelfolder,'SaveIntermediateResults');
            OCMATFINITCONT.basicresultfilename=fullocmatfile(modelfolder,'SaveIntermediateResultsGlobalVariable');
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4modelequilibrium.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4modelequilibrium.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        [numrow numcol]=size(xout);
        switch fieldname
            case {'modelequilibrium','modelequilibriumu2c'}
                ocC.y=xout(1:MODELINFO.MATCONTEQUILIBRIUM.nphase,1:numcol);
                ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
                ocC.arcposition=[1 numcol];
                ocC.modelname=modelname(ocObj);
                ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
                ocC.userinfo.varyparameterindex=MODELINFO.MATCONTEQUILIBRIUM.ActiveParams;
                ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
                if isfield(OCMATFINITCONT,'combinationfunction') && ~isempty(OCMATFINITCONT.combinationfunction)
                    ocC.userinfo.combinationfunction=OCMATFINITCONT.combinationfunction;
                end
            case 'modelequilibriummf'
                ocC.y([OCMATFINITCONT.equilibriummfcoord';OCMATFINITCONT.equilibriummffreecoord],:)=xout(1:MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
                ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
                ocC.arcposition=[1 numcol];
                ocC.modelname=modelname(ocObj);
                ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
                ocC.userinfo.varyparameterindex=MODELINFO.MATCONTEQUILIBRIUM.ActiveParams;
                ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTEQUILIBRIUM.nphase+1,1:numcol);
            case 'modelequilibriumc'
                ocC.y=xout(OCMATFINITCONT.depvarcoord,1:numcol);
                ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
                ocC.arcposition=[1 numcol];
                ocC.modelname=modelname(ocObj);
                ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
                ocC.userinfo.varyparameterindex=[OCMATFINITCONT.freeparameterindex;MODELINFO.MATCONTEQUILIBRIUM.ActiveParams];
                ocC.userinfo.varyparametervalue=xout(OCMATFINITCONT.freeparametercoord+[1 0],1:numcol);
        end
        ocC.userinfo.tangent=vout;
        ocC.userinfo.functioneval=fout;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationSolution=occurve(ocC);

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';

    case {'modellimitpoint'}
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            OCMATFINITCONTORIGINAL=[];
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4modellimitpoint.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4modellimitpoint.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        [numrow numcol]=size(xout);
        ocC.y=xout(1:MODELINFO.MATCONTEQUILIBRIUM.nphase,1:numcol);
        ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
        ocC.arcposition=[1 numcol];
        ocC.modelname=modelname(ocObj);
        ocC.modelparameter=MODELINFO.MATCONTEQUILIBRIUM.P0;
        ocC.userinfo.varyparameterindex=MODELINFO.MATCONTEQUILIBRIUM.ActiveParams;
        ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTEQUILIBRIUM.nphase+1:end,1:numcol);

        ocC.userinfo.tangent=vout;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationSolution=occurve(ocC);

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';
    case {'userfunction'}
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            OCMATFINITCONTORIGINAL=[];
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4userfunction.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4userfunction.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        [numrow numcol]=size(xout);
        ocC.y=xout;
        ocC.arcarg=MODELINFO.OCMATUSERFUNC.arcarg(ones(1,numcol));
        ocC.arcposition=[1 numcol];
        ocC.modelname=modelname(ocObj);
        ocC.modelparameter=MODELINFO.OCMATUSERFUNC.P0;

        ocC.userinfo.tangent=vout;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationSolution=occurve(ocC);

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';

    case 'modelhopf'
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            OCMATFINITCONTORIGINAL=[];
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4modelhopf.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4modelhopf.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        [numrow numcol]=size(xout);
        ocC.y=xout(1:MODELINFO.MATCONTHOPF.nphase,1:numcol);
        ocC.arcarg=OCMATFINITCONT.arcarg(ones(1,numcol));
        ocC.arcposition=[1 numcol];
        ocC.modelname=modelname(ocObj);
        ocC.modelparameter=MODELINFO.MATCONTHOPF.P0;
        ocC.userinfo.varyparameterindex=MODELINFO.MATCONTHOPF.ActiveParams;
        ocC.userinfo.varyparametervalue=xout(MODELINFO.MATCONTHOPF.nphase+(1:2),1:numcol);
        ocC.userinfo.kvalue=xout(MODELINFO.MATCONTHOPF.nphase+3,1:numcol);
        
        ocC.userinfo.tangent=vout;
        ocC.userinfo.functioneval=fout;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationSolution=occurve(ocC);
        
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';
        
    case 'modellimitcycle'
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            OCMATFINITCONTORIGINAL=[];
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4modellimitcycle.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4modellimitcycle.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        totalmesh=fout(1:MODELINFO.MATCONTLIMITCYCLE.ntst+1,:);
        if isfield(MODELINFO.MATCONTLIMITCYCLE,'multipliers') && ~isempty(MODELINFO.MATCONTLIMITCYCLE.multipliers)
            linearization=fout(MODELINFO.MATCONTLIMITCYCLE.ntst+MODELINFO.MATCONTLIMITCYCLE.nphase+2:end,:);
        else
            linearization=[];
        end
        ntstcolp1=MODELINFO.MATCONTLIMITCYCLE.ntstcol+1;
        nonemptyidx=find(sum(abs(xout)));
        modelnme=modelname(ocObj);
        counter=0;
        for ii=nonemptyidx
            counter=counter+1;
            trj(counter).y=reshape(xout(MODELINFO.MATCONTLIMITCYCLE.coords,ii),MODELINFO.MATCONTLIMITCYCLE.nphase,MODELINFO.MATCONTLIMITCYCLE.tps);
            msh=totalmesh(:,ii).';
            cumsum([msh(1:end-1);repmat(diff(msh)/MODELINFO.MATCONTLIMITCYCLE.ncol,MODELINFO.MATCONTLIMITCYCLE.ncol-1,1)],1);
            trj(counter).x=[ans(:).' 1];
            trj(counter).x0=0;
            trj(counter).timehorizon=xout(MODELINFO.MATCONTLIMITCYCLE.PeriodIdx,ii);
            trj(counter).arcarg=OCMATFINITCONT.arcarg;
            trj(counter).arcinterval=[0 xout(MODELINFO.MATCONTLIMITCYCLE.PeriodIdx,ii)];
            trj(counter).arcposition=[1;ntstcolp1];
            trj(counter).modelparameter=MODELINFO.MATCONTLIMITCYCLE.P0;
            trj(counter).modelparameter(MODELINFO.MATCONTLIMITCYCLE.ActiveParams)=xout(MODELINFO.MATCONTLIMITCYCLE.pars(end),ii);
            trj(counter).modelname=modelnme;
            trj(counter).solverinfo.tmesh=totalmesh(:,ii).';
            trj(counter).solverinfo.coeff=xout(:,ii);
            trj(counter).solverinfo.tangent=vout(:,ii);
            trj(counter).solverinfo.ntst=MODELINFO.MATCONTLIMITCYCLE.ntst;
            trj(counter).solverinfo.ncol=MODELINFO.MATCONTLIMITCYCLE.ncol;
            if ~isempty(linearization)
                trj(counter).linearization=reshape(linearization(:,ii),MODELINFO.MATCONTLIMITCYCLE.nphase,MODELINFO.MATCONTLIMITCYCLE.nphase);
            else
                trj(counter).linearization=[];
            end
        end
        ocResultStruct.ContinuationSolution=trj;
        if ~counter
            ocResultStruct=[];
            ocResultFieldName='';
            return
        end
        LastocObj=changeparametervalue(ocObj,trj(counter).modelparameter);
        perStruct.octrajectory=trj(counter);
        perStruct.period=trj(counter).arcinterval(2);

        ocResultStruct.LastLimitCycle=dynprimitive(perStruct,LastocObj);
        ocResultStruct.LastModel=LastocObj;
        ocResultStruct.ContinuationInformation=sout;

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';

    case 'modellimitpointcycle'
        global OCMATFINITCONT

        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATFINITCONT)
            OCMATFINITCONTORIGINAL=OCMATFINITCONT;
        else
            OCMATFINITCONTORIGINAL=[];
        end
        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATFINITCONT.basicresultfilename '4limitpointcycle.mat'];
        globalvarfile=[OCMATFINITCONT.basicglobalvarfilename '4limitpointcycle.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        totalmesh=fout(1:MODELINFO.MATCONTLIMITCYCLE.ntst+1,:);
        if isfield(MODELINFO.MATCONTLIMITCYCLE,'multipliers') && ~isempty(MODELINFO.MATCONTLIMITCYCLE.multipliers)
            linearization=fout(MODELINFO.MATCONTLIMITCYCLE.ntst+MODELINFO.MATCONTLIMITCYCLE.nphase+2:end,:);
        else
            linearization=[];
        end
        ntstcolp1=MODELINFO.MATCONTLIMITCYCLE.ntstcol+1;
        nonemptyidx=find(sum(abs(xout)));
        modelnme=modelname(ocObj);
        counter=0;
        for ii=nonemptyidx
            counter=counter+1;
            trj(counter).y=reshape(xout(MODELINFO.MATCONTLIMITCYCLE.coords,ii),MODELINFO.MATCONTLIMITCYCLE.nphase,MODELINFO.MATCONTLIMITCYCLE.tps);
            msh=totalmesh(:,ii).';
            cumsum([msh(1:end-1);repmat(diff(msh)/MODELINFO.MATCONTLIMITCYCLE.ncol,MODELINFO.MATCONTLIMITCYCLE.ncol-1,1)],1);
            trj(counter).x=[ans(:).' 1];
            trj(counter).x0=0;
            trj(counter).timehorizon=xout(MODELINFO.MATCONTLIMITCYCLE.PeriodIdx,ii);
            trj(counter).arcarg=OCMATFINITCONT.arcarg;
            trj(counter).arcinterval=[0 xout(MODELINFO.MATCONTLIMITCYCLE.PeriodIdx,ii)];
            trj(counter).arcposition=[1;ntstcolp1];
            trj(counter).modelparameter=MODELINFO.MATCONTLIMITCYCLE.P0;
            trj(counter).modelparameter(MODELINFO.MATCONTLIMITCYCLE.ActiveParams)=xout(MODELINFO.MATCONTLIMITCYCLE.pars(end-1:end),ii);
            trj(counter).modelname=modelnme;
            trj(counter).solverinfo.tmesh=totalmesh(:,ii).';
            trj(counter).solverinfo.coeff=xout(:,ii);
            trj(counter).solverinfo.tangent=vout(:,ii);
            trj(counter).solverinfo.ntst=MODELINFO.MATCONTLIMITCYCLE.ntst;
            trj(counter).solverinfo.ncol=MODELINFO.MATCONTLIMITCYCLE.ncol;
            if ~isempty(linearization)
                trj(counter).linearization=reshape(linearization(:,ii),MODELINFO.MATCONTLIMITCYCLE.nphase,MODELINFO.MATCONTLIMITCYCLE.nphase);
            else
                trj(counter).linearization=[];
            end
        end
        ocResultStruct.ContinuationSolution=trj;
        if ~counter
            ocResultStruct=[];
            ocResultFieldName='';
            return
        end
        LastocObj=changeparametervalue(ocObj,trj(counter).modelparameter);
        perStruct.octrajectory=trj(counter);
        perStruct.period=trj(counter).arcinterval(2);

        ocResultStruct.LastLimitCycle=dynprimitive(perStruct,LastocObj);
        ocResultStruct.LastModel=LastocObj;
        ocResultStruct.ContinuationInformation=sout;

        ocResultStruct.ContinuationClassification=fieldname;
        ocResultFieldName='MatContContinuation';

    case {'limitextremal','limitextremal2lc','limitextremal4ft'}
        global OCMATCONT OCMATLSC

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATLSCORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)
        mainfunch=feval(str2func(fieldname));
        % store value of the global variables OCMATCONT and OCMATLSC to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATLSC.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATLSC'' changed.\n')
        OCMATLSC=MODELINFO.OCMATLSC;

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                DataAdaptation(bvpout(ii).tmesh);
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                if any(strcmp(fieldname,{'limitextremal','limitextremal2lc'}))
                    limitpointStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                    limitpointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                    limitpointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                end
            end
        end
        limitpointStruct.arcposition=[];
        if any(strcmp(fieldname,{'limitextremal','limitextremal2lc'}))
            limitpointStruct.arcposition=[];
            ocEP=dynprimitive(OCMATLSC.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATLSC.linearization);
        end
        if ~isempty(soln)
            if any(strcmp(fieldname,{'limitextremal','limitextremal2lc'}))
                ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,ocObj),ocEP);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObj),ocEP);
            elseif any(strcmp(fieldname,{'limitextremal4ft'}))
                ocResultStruct.NonadmissibleSolution=octrajectory(soln,ocObj);
                ocResultStruct.ExtremalSolution=octrajectory(sol,ocObj);
            end
        else
            if ~isempty(soltarget)
                ocTrj=octrajectory(soltarget,ocObj);
            else
                ocTrj=octrajectory(ocExStruct(end),ocObj);
            end
            if any(strcmp(fieldname,{'limitextremal','limitextremal2lc'}))
                ocResultStruct.ExtremalSolution=ocasymptotic(ocTrj,ocEP);
            elseif any(strcmp(fieldname,{'limitextremal4ft'}))
                ocResultStruct.ExtremalSolution=ocTrj;
            end
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        if any(strcmp(fieldname,{'limitextremal','limitextremal2lc','limitextremal4ft'}))
            ocResultStruct.LimitPointCurve=occurve(limitpointStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if any(strcmp(fieldname,{'limitextremal','limitextremal2lc'}))
            ocResultStruct.LimitSet=ocEP;
        end
        if any(strcmp(fieldname,{'limitextremal4ft'})) && OCMATLSC.targettype==3
            ocObjtmp=changeparametervalue(ocObj,OCMATLSC.parameterindex,ocResultStruct.ExtremalSolution.solverinfo.parameters(end));
            ocResultStruct.LastModel=ocObjtmp;
        end
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATLSCORIGINAL,OCBVPORIGINAL);
        

    case {'limitextremalp'}
        global OCMATCONT OCMATLSC OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATLSCORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));


        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if isempty(soln) && isempty(soltarget)
            soltarget=ocExStruct(end);
        end
        if ~isempty(soln)
            ocObjtmp=changeparametervalue(ocObj,OCMATLSC.varyparameterindex,soln.solverinfo.parameters(end));
            equcoord=soln.solverinfo.equilibriumcoord;
            hatx.y=soln.solverinfo.parameters(equcoord);
            hatx.x=0;
            hatx.arcarg=OCMATCONT.HE.arcarg(end);
            hatx.linearization=jacobian(ocObjtmp,hatx);
            ocEP=dynprimitive(hatx);
            ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,ocObjtmp),ocEP);
            if isempty(sol)
                ocResultStruct.ExtremalSolution=ocasymptotic([]);
            else
                ocObjtmp=changeparametervalue(ocObj,OCMATLSC.varyparameterindex,soln.solverinfo.parameters(end));
                equcoord=soln.solverinfo.equilibriumcoord;
                hatx.y=soln.solverinfo.parameters(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                hatx.linearization=jacobian(ocObjtmp,hatx);
                ocEP=dynprimitive(hatx);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObjtmp),ocEP);
            end
        else
            ocObjtmp=changeparametervalue(ocObj,OCMATLSC.varyparameterindex,soltarget.solverinfo.parameters(end));
            equcoord=soltarget.solverinfo.equilibriumcoord;
            hatx.y=soltarget.solverinfo.parameters(equcoord);
            hatx.x=0;
            hatx.arcarg=OCMATCONT.HE.arcarg(end);
            hatx.linearization=jacobian(ocObjtmp,hatx);
            ocEP=dynprimitive(hatx);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,ocObjtmp),ocEP);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATLSCORIGINAL,OCBVPORIGINAL);

 
    case {'indifferencesolution4mm'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));


        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end

        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        %         if any(strcmp(fieldname,{'indifferencesolution4mm'})) && OCMATINDIF.targettype==3
        %             sol=ocExStruct(end);
        %             parnew=sol.modelparameter;
        %             mmObjtmp=changeparametervalue(mmObj,parnew);
        %             ocResultStruct.LastModel=mmObjtmp;
        %         end
        %ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);
       
    case {'indifferencesolution','indifferencesolution4emf','indifferencesolution4per','indifferencedistribution','indifferencesolution4ft','indifferencesolution4ae2ftae','indifferencesolution4ep_ft','indifferencesolutionep'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));


        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                indifferencepointStruct.y(OCMATINDIF.statecostatecoordinate,counter)=ocExStruct(ii).y(OCMATINDIF.statecostatecoordinate,1);
                indifferencepointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                indifferencepointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        indifferencepointStruct.arcposition=[];

        %         ocEP=dynprimitive(OCMATINDIF.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATINDIF.linearization);
        if ~isempty(soln)
            if any(strcmp(fieldname,{'indifferencesolution','indifferencedistribution','indifferencesolution4ae2ftae'}))
                switch soln.solver
                    case 'gbvp4c'
                        ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                        ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                    otherwise
                        ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,OCMATINDIF.saddlepoint,OCMATINDIF.linearization);
                        ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,OCMATINDIF.saddlepoint,OCMATINDIF.linearization);
                end
            elseif strcmp(fieldname,'indifferencesolution4per')
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
            elseif strcmp(fieldname,'indifferencesolution4emf')
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,soln.solverinfo.saddlepoint,OCMATINDIF.linearization);
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,sol.solverinfo.saddlepoint,OCMATINDIF.linearization);
            elseif strcmp(fieldname,'indifferencesolution4ft')
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln);
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol);
            elseif strcmp(fieldname,'indifferencesolutionep') 
                ocObjtmp=changeparametervalue(ocObj,OCMATINDIF.freeparameterindex,sol.solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                hatx.y=sol.solverinfo.parameters(sol.solverinfo.equilibriumcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                hatx.linearization=jacobian(ocObjtmp,hatx);
                ocEP=dynprimitive(hatx);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,ocObjtmp),ocEP);
                ocObjtmp=changeparametervalue(ocObj,OCMATINDIF.freeparameterindex,soln.solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                hatx.y=soln.solverinfo.parameters(sol.solverinfo.equilibriumcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(end);
                hatx.linearization=jacobian(ocObjtmp,hatx);
                ocEP=dynprimitive(hatx);
                ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,ocObjtmp),ocEP);
            elseif any(strcmp(fieldname,{'indifferencesolution4ep_ft','indifferencesolutionep'}))
                parnew=soln.modelparameter;
                parnew(OCMATINDIF.parameterindex)=soln.solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                ocObjtmp=changeparametervalue(ocObj,parnew);
                J=cell(1,OCMATINDIF.indifferenceorder);
                saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                for ii=1:OCMATINDIF.indifferenceorder
                    if OCMATINDIF.ocasymptotic(ii)
                        equcoord=soln.solverinfo.equilibriumcoord{ii};
                    else
                        equcoord=[];
                    end
                    if ~isempty(equcoord)
                        hatx.y=soln.solverinfo.parameters(equcoord);
                        hatx.x=0;
                        hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                        J{ii}=jacobian(ocObjtmp,hatx);
                        saddlepoint{ii}=hatx.y;
                    end
                end
                ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATINDIF,saddlepoint,J);
                parnew=sol.modelparameter;
                parnew(OCMATINDIF.parameterindex)=sol.solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                ocObjtmp=changeparametervalue(ocObj,parnew);
                J=cell(1,OCMATINDIF.indifferenceorder);
                saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                for ii=1:OCMATINDIF.indifferenceorder
                    if OCMATINDIF.ocasymptotic(ii)
                        equcoord=sol.solverinfo.equilibriumcoord{ii};
                    else
                        equcoord=[];
                    end
                    if ~isempty(equcoord)
                        hatx.y=sol.solverinfo.parameters(equcoord);
                        hatx.x=0;
                        hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                        J{ii}=jacobian(ocObjtmp,hatx);
                        saddlepoint{ii}=hatx.y;
                    end
                end
                ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATINDIF,saddlepoint,J);
            end
        else
            if ~isempty(soltarget)
                if any(strcmp(fieldname,{'indifferencesolution','indifferencedistribution','indifferencesolution4ae2ftae'}))
                    switch soltarget.solver
                        case 'gbvp4c'
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                        otherwise
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,OCMATINDIF.saddlepoint,OCMATINDIF.linearization);
                    end
                elseif strcmp(fieldname,'indifferencesolution4per')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4emf')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,soltarget.solverinfo.saddlepoint,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4ft')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget);
                elseif strcmp(fieldname,'indifferencesolutionep')
                    ocObjtmp=changeparametervalue(ocObj,OCMATINDIF.freeparameterindex,soltarget.solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                    hatx.y=soltarget.solverinfo.parameters(soltarget.solverinfo.equilibriumcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(end);
                    hatx.linearization=jacobian(ocObjtmp,hatx);
                    ocEP=dynprimitive(hatx);
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,ocObjtmp),ocEP);
                elseif strcmp(fieldname,'indifferencesolution4ep_ft')
                    parnew=soltarget.modelparameter;
                    parnew(OCMATINDIF.parameterindex)=soltarget.solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                    ocObjtmp=changeparametervalue(ocObj,parnew);
                    J=cell(1,OCMATINDIF.indifferenceorder);
                    saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                    for ii=1:OCMATINDIF.indifferenceorder
                        if OCMATINDIF.ocasymptotic(ii)
                            equcoord=soltarget.solverinfo.equilibriumcoord{ii};
                        else
                            equcoord=[];
                        end
                        if ~isempty(equcoord)
                            hatx.y=soltarget.solverinfo.parameters(equcoord);
                            hatx.x=0;
                            hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                            J{ii}=jacobian(ocObjtmp,hatx);
                            saddlepoint{ii}=hatx.y;
                        end
                    end
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATINDIF,saddlepoint,J);
                end
            else
                if any(strcmp(fieldname,{'indifferencesolution','indifferencedistribution','indifferencesolution4ae2ftae'}))
                    switch ocExStruct(end).solver
                        case 'gbvp4c'
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                        otherwise
                            ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                    end
                elseif strcmp(fieldname,'indifferencesolution4per')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,OCMATINDIF.limitset,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4emf')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,ocExStruct(end).solverinfo.saddlepoint,OCMATINDIF.linearization);
                elseif strcmp(fieldname,'indifferencesolution4ft')
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end));
                elseif strcmp(fieldname,'indifferencesolutionep')
                    ocObjtmp=changeparametervalue(ocObj,OCMATINDIF.freeparameterindex,ocExStruct(end).solverinfo.parameters(OCMATINDIF.freeparametercoordinate));
                    hatx.y=ocExStruct(end).solverinfo.parameters(ocExStruct(end).solverinfo.equilibriumcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(end);
                    hatx.linearization=jacobian(ocObjtmp,hatx);
                    ocEP=dynprimitive(hatx);
                    ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),ocObjtmp),ocEP);
                elseif strcmp(fieldname,'indifferencesolution4ep_ft')
                    parnew=ocExStruct(end).modelparameter;
                    parnew(OCMATINDIF.parameterindex)=ocExStruct(end).solverinfo.parameters(OCMATINDIF.parametervaluecoord);
                    ocObjtmp=changeparametervalue(ocObj,parnew);
                    J=cell(1,OCMATINDIF.indifferenceorder);
                    saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
                    for ii=1:OCMATINDIF.indifferenceorder
                        if OCMATINDIF.ocasymptotic(ii)
                            equcoord=ocExStruct(end).solverinfo.equilibriumcoord{ii};
                        else
                            equcoord=[];
                        end
                        if ~isempty(equcoord)
                            hatx.y=ocExStruct(end).solverinfo.parameters(equcoord);
                            hatx.x=0;
                            hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                            J{ii}=jacobian(ocObjtmp,hatx);
                            saddlepoint{ii}=hatx.y;
                        end
                    end
                    ocResultStruct.ExtremalSolution=sol2ocmultipath(ocExStruct(end),OCMATINDIF,saddlepoint,J);
                end
                if strcmp(fieldname,'indifferencesolution4ft')
                    ocResultStruct.NonadmissibleSolution=octrajectory([]);
                else
                    ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
                end
            end
        end
        if any(strcmp(fieldname,{'indifferencesolution','indifferencesolution4emf','indifferencesolution4per','indifferencedistribution','indifferencesolution4ft','indifferencesolution4ae2ftae','indifferencesolution4ep_ft','indifferencesolutionep'}))
            ocResultStruct.IndifferencePointCurve=occurve(indifferencepointStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if any(strcmp(fieldname,{'indifferencesolution4ft'})) && OCMATINDIF.continuationtype==3
            ocObjtmp=changeparametervalue(ocObj,OCMATINDIF.parameterindex,ocExStruct(end).solverinfo.parameters(OCMATINDIF.parametercoord));
            ocResultStruct.LastModel=ocObjtmp;
        elseif any(strcmp(fieldname,{'indifferencesolution4ep_ft'}))
            ocResultStruct.LastModel=ocObjtmp;
        end
        %ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);


    case {'indifferencesolutionp','indifferencesolution2lc'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));


        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if isempty(soln) && isempty(soltarget)
            soltarget=ocExStruct(end);
        end
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=indifferencesolution2ocmultipath(ocObj,soln);
            ocResultStruct.ExtremalSolution=indifferencesolution2ocmultipath(ocObj,sol);
        else
            ocResultStruct.ExtremalSolution=indifferencesolution2ocmultipath(ocObj,soltarget);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);

    case {'heteroclinic','heterocliniccyc','heteroclinicep2ft','homoclinic'}
        global OCMATCONT OCMATHET OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATHETORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        %OCMATHET.hetorder=2;
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if isempty(soltarget)
            soltarget=ocExStruct(end);
        end
        if 0%~isempty(soln)
            ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATHET,OCMATHET.saddlepoint,OCMATHET.linearization);
            ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATHET,OCMATHET.saddlepoint,OCMATHET.linearization);
        else
            parnew=soltarget.modelparameter;
            parnew(OCMATHET.parameterindex)=soltarget.solverinfo.parameters(OCMATHET.parametervaluecoord);
            ocObjtmp=changeparametervalue(ocObj,parnew);
            J=cell(1,OCMATHET.hetorder);
            saddlepoint=cell(1,OCMATHET.hetorder);
            for ii=1:OCMATHET.hetorder
                if any(strcmp(fieldname,{'heteroclinic'}))
                    equcoord=soltarget.solverinfo.equilibriumcoord{ii};
                elseif any(strcmp(fieldname,{'heteroclinicep2ft'}))
                    if OCMATHET.ocasymptotic(ii)
                        equcoord=soltarget.solverinfo.equilibriumcoord{ii};
                    else
                        equcoord=[];
                    end
                elseif any(strcmp(fieldname,{'heterocliniccyc'}))
                    equcoord=soltarget.solverinfo.equilibriumcoord{OCMATHET.equilibriumidx(ii)};
                else
                    equcoord=soltarget.solverinfo.equilibriumcoord;
                end
                if ~isempty(equcoord)
                    hatx.y=soltarget.solverinfo.parameters(equcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(OCMATHET.arccoord{ii}(end));
                    J{ii}=jacobian(ocObjtmp,hatx);
                    saddlepoint{ii}=hatx.y;
                end
            end
            OCMATHET.indifferenceorder=OCMATHET.hetorder;
            if any(strcmp(fieldname,{'heteroclinic','heterocliniccyc','heteroclinicep2ft'}))
                ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATHET,saddlepoint,J);
            end
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATHETORIGINAL,OCBVPORIGINAL);


    case {'extremalmp2ep'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        OCMATAE.hetorder=2;
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if isempty(soltarget)
            soltarget=ocExStruct(end);
        end
        if 0%~isempty(soln)
            ocResultStruct.NonadmissibleSolution=sol2ocmultipath(soln,OCMATAE,OCMATAE.saddlepoint,OCMATAE.linearization);
            ocResultStruct.ExtremalSolution=sol2ocmultipath(sol,OCMATAE,OCMATAE.saddlepoint,OCMATAE.linearization);
        else
            parnew=soltarget.modelparameter;
            parnew(OCMATAE.parameterindex)=soltarget.solverinfo.parameters(OCMATAE.parametervaluecoord);
            ocObjtmp=changeparametervalue(ocObj,parnew);
            J=cell(1,OCMATAE.multorder);
            saddlepoint=cell(1,OCMATAE.multorder);
            for ii=1:OCMATAE.multorder
                calcequilibrium=1;
                equcoord=soltarget.solverinfo.equilibriumcoord;
                if calcequilibrium
                    hatx.y=soltarget.solverinfo.parameters(equcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(OCMATAE.arccoord{ii}(end));
                    J{ii}=jacobian(ocObjtmp,hatx);
                    saddlepoint{ii}=hatx.y;
                end
            end
            OCMATAE.indifferenceorder=OCMATAE.multorder;
            ocResultStruct.ExtremalSolution=sol2ocmultipath(soltarget,OCMATAE,saddlepoint,J);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'limitcycle','limitcycleuser','limitcycleuserI','limitlimitcycle'}
        global OCMATCONT OCMATLC OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        OCBVP.F=[];
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            catch
                lasterr
            end
        end
        if ~isempty(soln)
            if any(strcmp(fieldname,{'limitcycle','limitcycleuser','limitlimitcycle'}))
                ocObjtmp=changeparametervalue(ocObj,soln.octrajectory.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            ocResultStruct.NonadmissibleLimitCycle=dynprimitive(soln,ocObjtmp);
            ocResultStruct.LimitCycle=dynprimitive(sol,ocObjtmp);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.LimitCycleMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                if any(strcmp(fieldname,{'limitcycle','limitcycleuser','limitlimitcycle'}))
                    ocObjtmp=changeparametervalue(ocObj,soltarget.octrajectory.modelparameter);
                else
                    ocObjtmp=ocObj;
                end
                ocResultStruct.LimitCycle=dynprimitive(soltarget,ocObjtmp);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.LimitCycleMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                if any(strcmp(fieldname,{'limitcycle','limitcycleuser','limitlimitcycle'}))
                    ocObjtmp=changeparametervalue(ocObj,ocExStruct(end).octrajectory.modelparameter);
                else
                    ocObjtmp=ocObj;
                end
                ocResultStruct.LimitCycle=dynprimitive(sollast);
                trajectory=ocExStruct(end).octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.LimitCycleMatlabFormat=matlabsol;
                end
            end
            ocResultStruct.NonadmissibleLimitCycle=dynprimitive([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'modelperiodicsol'}
        global OCMATCONT OCMATPS OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if ~isempty(soln)
            ocObjtmp=changeparametervalue(ocObj,OCMATPS.varyparameterindex,soln.octrajectory.solverinfo.continuationparameter);
            ocResultStruct.NonadmissiblePeriodicSolution=dynprimitive(soln,ocObjtmp);
            ocResultStruct.PeriodicSolution=dynprimitive(sol,ocObjtmp);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocObjtmp=changeparametervalue(ocObj,OCMATPS.varyparameterindex,soltarget.octrajectory.solverinfo.continuationparameter);
                ocResultStruct.PeriodicSolution=dynprimitive(soltarget,ocObjtmp);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
                end
            else
                %                 DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                ocObjtmp=changeparametervalue(ocObj,OCMATPS.varyparameterindex,ocExStruct(end).octrajectory.solverinfo.continuationparameter);
                %                 ocResultStruct.PeriodicSolution=dynprimitive(soltarget,ocObjtmp);
                %                 trajectory=ocExStruct(end).octrajectory;
                %                 matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                %                 if ~isempty(matlabsol)
                %                     ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
                %                 end
            end
            ocResultStruct.NonadmissiblePeriodicSolution=dynprimitive([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=ocObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'periodicsol_T'}
        global OCMATCONT OCMATPS OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if ~isempty(soln)
            ocResultStruct.NonadmissiblePeriodicSolution=dynprimitive(soln,ocObj);
            ocResultStruct.PeriodicSolution=dynprimitive(sol,ocObj);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocResultStruct.PeriodicSolution=dynprimitive(soltarget,ocObj);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
                end
            else
                %                 DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                %                 ocResultStruct.PeriodicSolution=dynprimitive(soltarget,ocObjtmp);
                %                 trajectory=ocExStruct(end).octrajectory;
                %                 matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                %                 if ~isempty(matlabsol)
                %                     ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
                %                 end
            end
            ocResultStruct.NonadmissiblePeriodicSolution=dynprimitive([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
    case {'extremal2per','extremal2lc'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(1:OCBVP.numode,counter)=ocExStruct(ii).y(1:OCBVP.numode,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        slMfStruct.arcposition=[];
        if ~isempty(soln)
                switch OCMATCONT.continuationtype
                    case 'parameter'
                        lcidx=find(sollast.solverinfo.solutionindex==2);
                        freepar=sollast.solverinfo.parameters;
                        arcposition=sollast.arcposition(:,lcidx);
                        ocTrj.x=sollast.x(arcposition(1,1):arcposition(2,end));
                        ocTrj.x=ocTrj.x-ocTrj.x(1);
                        ocTrj.y=sollast.y(:,arcposition(1,1):arcposition(2,end));
                        ocTrj.arcposition=arcposition-arcposition(1,1)+1;
                        ocTrj.arcarg=OCMATAE.LC.arcarg;
                        ocTrj.arcinterval=[0 freepar([OCMATAE.LC.switchtimecoordinate(:).' OCMATAE.LC.periodcoordinate]).'];
                        ocTrj.modelparameter=sollast.modelparameter;
                        ocTrj.solverinfo=sollast.solverinfo;
                        ocLC.octrajectory=ocTrj;
                        ocLC.period=ocTrj.arcinterval(end);
                        ocLC.octrajectory.linearization=sollast.solverinfo.monodromy;
                        ocLC=dynprimitive(ocLC);

                        lcidx=find(soln.solverinfo.solutionindex==2);
                        freepar=soln.solverinfo.parameters;
                        arcposition=soln.arcposition(:,lcidx);
                        ocTrj.x=soln.x(arcposition(1,1):arcposition(2,end));
                        ocTrj.x=ocTrj.x-ocTrj.x(1);
                        ocTrj.y=soln.y(:,arcposition(1,1):arcposition(2,end));
                        ocTrj.arcposition=arcposition-arcposition(1,1)+1;
                        ocTrj.arcarg=OCMATAE.LC.arcarg;
                        ocTrj.arcinterval=[0 freepar([OCMATAE.LC.switchtimecoordinate(:).' OCMATAE.LC.periodcoordinate]).'];
                        ocTrj.modelparameter=soln.modelparameter;
                        ocTrj.solverinfo=soln.solverinfo;
                        ocLCN.octrajectory=ocTrj;
                        ocLCN.period=ocTrj.arcinterval(end);
                        ocLCN.octrajectory.linearization=soln.solverinfo.monodromy;
                        ocLCN=dynprimitive(ocLCN);

                        
                        lcidx=find(soln.solverinfo.solutionindex==1);
                        freepar=soln.solverinfo.parameters;
                        arcposition=soln.arcposition(:,lcidx);
                        ocTrjN.x=soln.x(arcposition(1,1):arcposition(2,end));
                        ocTrjN.y=soln.y(:,arcposition(1,1):arcposition(2,end));
                        ocTrjN.arcposition=arcposition-arcposition(1,1)+1;
                        ocTrjN.arcarg=OCMATAE.TRJ.arcarg;
                        ocTrjN.arcinterval=[0 freepar([OCMATAE.TRJ.switchtimecoordinate(:).' OCMATAE.TRJ.endtimecoordinate]).'];
                        ocTrjN.modelparameter=soln.modelparameter;
                        ocTrjN.solverinfo=soln.solverinfo;
  
                        lcidx=find(sollast.solverinfo.solutionindex==1);
                        freepar=sollast.solverinfo.parameters;
                        arcposition=sollast.arcposition(:,lcidx);
                        ocTrj.x=sollast.x(arcposition(1,1):arcposition(2,end));
                        ocTrj.y=sollast.y(:,arcposition(1,1):arcposition(2,end));
                        ocTrj.arcposition=arcposition-arcposition(1,1)+1;
                        ocTrj.arcarg=OCMATAE.TRJ.arcarg;
                        ocTrj.arcinterval=[0 freepar([OCMATAE.TRJ.switchtimecoordinate(:).' OCMATAE.TRJ.endtimecoordinate]).'];
                        ocTrj.modelparameter=sollast.modelparameter;
                        ocTrj.solverinfo=sollast.solverinfo;

                        ocObj=changeparametervalue(ocObj,ocTrj.modelparameter);
                        ocObjN=changeparametervalue(ocObj,ocTrjN.modelparameter);
                        
                    otherwise
                        ocLC=OCMATAE.limitset;
                        ocLCN=OCMATAE.limitset;
                        ocTrj=sollast;
                        ocTrjN=soln;
                        ocObjN=ocObj;
                end
            ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(ocTrjN,ocObjN),ocLCN);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocTrj,ocObj),ocLC);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            if isempty(soltarget)
                soltarget=sollast;
            end
            if ~isempty(soltarget)
                switch OCMATCONT.continuationtype
                    case 'parameter'
                        lcidx=find(soltarget.solverinfo.solutionindex==2);
                        freepar=soltarget.solverinfo.parameters;
                        arcposition=soltarget.arcposition(:,lcidx);
                        ocTrj.x=soltarget.x(arcposition(1,1):arcposition(2,end));
                        ocTrj.x=ocTrj.x-ocTrj.x(1);
                        ocTrj.y=soltarget.y(:,arcposition(1,1):arcposition(2,end));
                        ocTrj.arcposition=arcposition-arcposition(1,1)+1;
                        ocTrj.arcarg=OCMATAE.LC.arcarg;
                        ocTrj.arcinterval=[0 freepar([OCMATAE.LC.switchtimecoordinate(:).' OCMATAE.LC.periodcoordinate]).'];
                        ocTrj.modelparameter=soltarget.modelparameter;
                        ocTrj.solverinfo=soltarget.solverinfo;
                        ocLC.octrajectory=ocTrj;
                        ocLC.period=ocTrj.arcinterval(end);
                        ocLC.octrajectory.linearization=soltarget.solverinfo.monodromy;
                        ocLC=dynprimitive(ocLC);
 
                        lcidx=find(soltarget.solverinfo.solutionindex==1);
                        freepar=soltarget.solverinfo.parameters;
                        arcposition=soltarget.arcposition(:,lcidx);
                        ocTrj.x=soltarget.x(arcposition(1,1):arcposition(2,end));
                        ocTrj.y=soltarget.y(:,arcposition(1,1):arcposition(2,end));
                        ocTrj.arcposition=arcposition-arcposition(1,1)+1;
                        ocTrj.arcarg=OCMATAE.TRJ.arcarg;
                        ocTrj.arcinterval=[0 freepar([OCMATAE.TRJ.switchtimecoordinate(:).' OCMATAE.TRJ.endtimecoordinate]).'];
                        ocTrj.modelparameter=soltarget.modelparameter;
                        ocTrj.solverinfo=soltarget.solverinfo;

                        ocObj=changeparametervalue(ocObj,ocTrj.modelparameter);
                    otherwise
                        ocLC=OCMATAE.limitset;
                        ocTrj=soltarget;
                end
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocTrj,ocObj),ocLC);
                matlabsol=transform2nativematlab(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).solverinfo.tmesh);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),ocObj),OCMATAE.limitset);
                matlabsol=transform2nativematlab(ocExStruct(end).solverinfo.tmesh,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            end
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        if any( strcmp(fieldname,{'extremal2per','extremal2lc'}))
            ocResultStruct.SliceManifold=occurve(slMfStruct);
        end
        
        
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocLC;
        if strcmp(OCMATCONT.continuationtype,'parameter')
            ocResultStruct.LastModel=ocObj;
        end
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'odesolution'}
        global OCMATCONT OCMATODEBVP OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=octrajectory(soln);
            ocResultStruct.ExtremalSolution=octrajectory(sol);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocResultStruct.ExtremalSolution=octrajectory(soltarget);
                matlabsol=transform2nativematlab(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).solverinfo.tmesh);
                ocResultStruct.ExtremalSolution=octrajectory(ocExStruct(end));
                matlabsol=transform2nativematlab(ocExStruct(end).solverinfo.tmesh,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            end
            ocResultStruct.NonadmissibleSolution=octrajectory([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationUserFunction=func2str(OCMATCONT.modelfunc);
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
    case {'extremal2ftae'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(counter).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        slMfStruct.arcposition=[];

        ocResultStruct.SliceManifold=occurve(slMfStruct);
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=octrajectory(soln);
            ocResultStruct.ExtremalSolution=octrajectory(sol);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocResultStruct.ExtremalSolution=octrajectory(soltarget);
                matlabsol=transform2nativematlab(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).solverinfo.tmesh);
                ocResultStruct.ExtremalSolution=octrajectory(ocExStruct(end));
                matlabsol=transform2nativematlab(ocExStruct(end).solverinfo.tmesh,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            end
            ocResultStruct.NonadmissibleSolution=octrajectory([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.ContinuationUserFunction=func2str(OCMATCONT.modelfunc);
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
        
%%%%%
%% solutions calculated by GRADIENT method
%%%%%

    case {'extremalgrad4ft','indifferencegradsolution','extremalgrad2ep'}
        global OCGRADCONT OCGRADSOL
        [globalvarfile,resultfile,OCGRADCONTORIGINAL,OCGRADSOLORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        contclass=OCGRADCONT.problem_func;
        conttype=OCGRADCONT.conttype;
        for ii=1:numel(gradout)
            try
                counter=counter+1;
                ocExStruct(counter)=mainfunch{22}(gradout(ii).coeff,gradout(ii).extremal,gradout(ii).tangent);
                ocExStruct(counter).extremal.gradientinfo=gradout(ii).gradientinfo;
            end
        end
        %OCMATCONT.meshadaptationflag=meshadaptationflag;
        if ~isempty(soltarget)
            sol=soltarget;
            ocResultStruct.ExtremalSolution=ocgradtrajectory(sol.extremal(1));
            switch contclass
                case 'extremalgrad2ep'
                    if ~strcmp(conttype,'parameter')
                        ocResultStruct.ExtremalSolution=ocgradasymptotic(ocResultStruct.ExtremalSolution,OCGRADSOL.equilibrium);
                    else
                        ocmatmsg('Not implemented yet')
                    end
                case 'indifferencegradsolution'
                    if strcmp(OCGRADCONT.trajectoryclass{1},'inf')
                        ocResultStruct.ExtremalSolution=ocgradasymptotic(ocResultStruct.ExtremalSolution,OCGRADSOL.equilibrium{1});
                    end
                    for ii=2:length(sol.extremal)
                        ocResultStruct.ExtremalSolution(ii)=ocgradtrajectory(sol.extremal(ii));
                        if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
                            ocResultStruct.ExtremalSolution(ii)=ocgradasymptotic(ocResultStruct.ExtremalSolution(ii),OCGRADSOL.equilibrium{ii});
                        end
                    end
            end
            ocResultStruct.NonadmissibleSolution=ocgradtrajectory([]);
        elseif ~isempty(soln)
            ocResultStruct.ExtremalSolution=ocgradtrajectory(sol.extremal(1));
            switch contclass
                case 'extremalgrad2ep'
                    if ~strcmp(conttype,'parameter')
                        ocResultStruct.ExtremalSolution=ocgradasymptotic(ocResultStruct.ExtremalSolution,OCGRADSOL.equilibrium);
                    else
                        ocmatmsg('Not implemented yet')
                    end
                case 'indifferencegradsolution'
                    if strcmp(OCGRADCONT.trajectoryclass{1},'inf')
                        ocResultStruct.ExtremalSolution=ocgradasymptotic(ocResultStruct.ExtremalSolution,OCGRADSOL.equilibrium{1});
                    end
                    for ii=2:length(sol.extremal)
                        ocResultStruct.ExtremalSolution(ii)=ocgradtrajectory(sol.extremal(ii));
                        if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
                            ocResultStruct.ExtremalSolution(ii)=ocgradasymptotic(ocResultStruct.ExtremalSolution(ii),OCGRADSOL.equilibrium{ii});
                        end
                    end
            end
            ocResultStruct.NonadmissibleSolution=ocgradtrajectory(soln.extremal(1));
            switch contclass
                case 'extremalgrad2ep'
                    if ~strcmp(conttype,'parameter')
                        ocResultStruct.NonadmissibleSolution=ocgradasymptotic(ocResultStruct.NonadmissibleSolution,OCGRADSOL.equilibrium);
                    else
                        ocmatmsg('Not implemented yet')
                    end
                case 'indifferencegradsolution'
                    if strcmp(OCGRADCONT.trajectoryclass{1},'inf')
                        ocResultStruct.NonadmissibleSolution=ocgradasymptotic(ocResultStruct.NonadmissibleSolution,OCGRADSOL.equilibrium{1});
                    end
                    for ii=2:length(sol.extremal)
                        ocResultStruct.NonadmissibleSolution(ii)=ocgradtrajectory(soln.extremal(ii));
                        if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
                            ocResultStruct.NonadmissibleSolution(ii)=ocgradasymptotic(ocResultStruct.NonadmissibleSolution(ii),OCGRADSOL.equilibrium{ii});
                        end
                    end
            end
        else
            sol=ocExStruct(end);
            ocResultStruct.ExtremalSolution=ocgradtrajectory(sol.extremal(1));
            switch contclass
                case 'extremalgrad2ep'
                    if ~strcmp(conttype,'parameter')
                        ocResultStruct.ExtremalSolution=ocgradasymptotic(ocResultStruct.ExtremalSolution,OCGRADSOL.equilibrium);
                    else
                        ocmatmsg('Not implemented yet')
                    end
                case 'indifferencegradsolution'
                    if strcmp(OCGRADCONT.trajectoryclass{1},'inf')
                        ocResultStruct.ExtremalSolution=ocgradasymptotic(ocResultStruct.ExtremalSolution,OCGRADSOL.equilibrium{1});
                    end
                    for ii=2:length(sol.extremal)
                        if strcmp(OCGRADCONT.trajectoryclass{ii},'inf')
                            ocResultStruct.ExtremalSolution(ii)=ocgradasymptotic(ocgradtrajectory(sol.extremal(ii)),OCGRADSOL.equilibrium{ii});
                        else
                            ocResultStruct.ExtremalSolution(ii)=ocgradtrajectory(sol.extremal(ii));
                        end
                    end
            end
            ocResultStruct.NonadmissibleSolution=ocgradtrajectory([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.GlobalInformation=OCGRADCONT;
        ocResultStruct.ModelInformation=OCGRADSOL;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCGRADCONTORIGINAL,OCGRADSOLORIGINAL);


    case {'extremaldae4ft'}
        global OCMATCONT OCMATFTE 
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATFTEORIGINAL]=loadglobal(ocObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);

        counter=0;
        %meshadaptationflag=OCMATCONT.meshadaptationflag;
        OCMATCONT.meshadaptationflag=1;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                %DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tcolmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                ocExStruct(counter).solverinfo.stepwidth=bvpout(ii).stepwidth;
                slMfStruct.y(:,counter)=ocExStruct(ii).y(:,1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                if OCMATFTE.objectivevaluecalc
                    slMfStruct.userinfo.objectivevalue(counter)=ocExStruct(ii).y(OCMATFTE.objectivevaluecoord,end);
                end
            end
        end
        slMfStruct.arcarg=[];
        slMfStruct.arcposition=[];
        %OCMATCONT.meshadaptationflag=meshadaptationflag;
        if ~isempty(soltarget)
            sol=soltarget;
            if strcmp(sol.solverinfo.continuationtype,'parameter')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            ocResultStruct.ExtremalSolution=octrajectory2ocdae(octrajectory(sol,ocObjtmp));
            ocResultStruct.NonadmissibleSolution=octrajectory([]);
        elseif ~isempty(soln)
            if strcmp(soln.solverinfo.continuationtype,'parameter')
                ocObjtmp=changeparametervalue(ocObj,soln.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            ocResultStruct.NonadmissibleSolution=octrajectory2ocdae(octrajectory(sol,ocObjtmp));
            if strcmp(sol.solverinfo.continuationtype,'parameter')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            ocResultStruct.ExtremalSolution=octrajectory2ocdae(octrajectory(sol,ocObjtmp));
        else
            sol=ocExStruct(end);
            if strcmp(sol.solverinfo.continuationtype,'parameter')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            else
                ocObjtmp=ocObj;
            end
            ocResultStruct.ExtremalSolution=octrajectory2ocdae(octrajectory(sol,ocObjtmp));
            ocResultStruct.NonadmissibleSolution=octrajectory([]);
        end
        ocResultStruct.SliceManifold=occurve(slMfStruct);
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if strcmp(sol.solverinfo.continuationtype,'parameter')
            ocResultStruct.ContinuationParameter=OCMATFTE.continuationindex;
            ocResultStruct.LastModel=ocObjtmp;
        end
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATFTEORIGINAL);

        
    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(ocObj,fieldname)
switch fieldname
    case {'extremal2ep','extremalt2ep','extremalp2ep','extremalc2ep','extremalp2epuser','extremale2ep','extremal2per','extremal2lc','extremal2emf','extremalt2emf','optisocline','extremal2inf','extremalt2inf','extremalp2inf','extremal2ftae','saddlepath2ep'}
        global OCMATCONT OCMATAE OCBVP
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATAE)
            OCMATAEORIGINAL=OCMATAE;
        else
            OCMATAEORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATAE)
            resultfile=[OCMATAE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATAE.basicglobalvarfilename '4' fieldname '.mat'];
        end

        if isempty(OCMATAE) || ~exist(resultfile,'file')
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
        end
        if isempty(OCMATAE) || ~exist(globalvarfile,'file')
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end
        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATAE'' changed.\n')
        OCMATAE=MODELINFO.OCMATAE;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATAEORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'extremal4ft','extremalt4ft','extremalp4ft','extremaltF4ft','extremalpF4ft','extremalF4ft','extremalt4inft'}
        global OCMATCONT OCBVP OCMATFTE
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end
        if ~isempty(OCMATFTE)
            OCMATFTEORIGINAL=OCMATFTE;
        else
            OCMATFTEORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATFTE)
            resultfile=[OCMATFTE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATFTE.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATFTE'' changed.\n')
        OCMATFTE=MODELINFO.OCMATFTE;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATFTEORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'limitextremal','limitextremalp','limitextremal2lc','limitextremal4ft'}
        global OCMATCONT OCMATLSC OCBVP

        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATLSC)
            OCMATLSCORIGINAL=OCMATLSC;
        else
            OCMATLSCORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end
        % mainfunch{22} function handle to format solution returned by the
        % bvp continuation process into structure that is consistent with
        % octrajectory class
        if any(strcmp(fieldname,{'limitextremal','limitextremal2lc'}))
            funch=feval(modelspecificfunc(ocObj,'4SaddlePathContinuation'));
        else
            funch=feval(modelspecificfunc(ocObj,'4FiniteHorizonPathContinuation'));
        end
        % funch{21} returns the basic model specific file names where
        % continuation information is stored

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        resultfile=[OCMATLSC.basicresultfilename '4' fieldname '.mat'];
        globalvarfile=[OCMATLSC.basicglobalvarfilename '4' fieldname '.mat'];

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end

        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATLSC'' changed.\n')
        OCMATLSC=MODELINFO.OCMATLSC;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATLSCORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'indifferencesolution','indifferencesolution4emf','indifferencesolution4per','indifferencedistribution','indifferencesolution4ft','indifferencesolutionp','indifferencesolution4ae2ftae','indifferencesolution4ep_ft','indifferencesolutionep','indifferencesolution4mm', ...
            'indifferencesolution4lcae','indifferencesolution2lc'}
        global OCMATCONT OCBVP OCMATINDIF
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end
        if ~isempty(OCMATINDIF)
            OCMATINDIFORIGINAL=OCMATINDIF;
        else
            OCMATINDIFORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATINDIF)
            resultfile=[OCMATINDIF.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATINDIF.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATINDIF'' changed.\n')
        OCMATINDIF=MODELINFO.OCMATINDIF;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATINDIFORIGINAL;
        varargout{3}=OCBVPORIGINAL;


    case {'heteroclinic','heterocliniccyc','heteroclinicep2ft','homoclinic','heteroclinicep2lc'}
        global OCMATCONT OCBVP OCMATHET
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end
        if ~isempty(OCMATHET)
            OCMATHETORIGINAL=OCMATHET;
        else
            OCMATHETORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATHET)
            resultfile=[OCMATHET.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATHET.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATHET'' changed.\n')
        OCMATHET=MODELINFO.OCMATHET;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATHETORIGINAL;
        varargout{3}=OCBVPORIGINAL;


    case {'extremalmp2ep'}
        global OCMATCONT OCBVP OCMATAE
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end
        if ~isempty(OCMATAE)
            OCMATAEORIGINAL=OCMATAE;
        else
            OCMATAEORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATAE)
            resultfile=[OCMATAE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATAE.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATAE'' changed.\n')
        OCMATAE=MODELINFO.OCMATAE;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATAEORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'limitcycle','limitcycleuser','limitcycleuserI','limitlimitcycle'}
        global OCMATCONT OCMATLC OCBVP
        % store value of the global variables OCMATCONT and OCMATLC to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATLC)
            OCMATLCORIGINAL=OCMATLC;
        else
            OCMATLCORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATLC)
            resultfile=[OCMATLC.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATLC.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATLC.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATLC'' changed.\n')
        OCMATLC=MODELINFO.OCMATLC;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATLCORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'modelperiodicsol','periodicsol_T'}
        global OCMATCONT OCMATPS OCBVP
        % store value of the global variables OCMATCONT and OCMATPS to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATPS)
            OCMATPSORIGINAL=OCMATPS;
        else
            OCMATPSORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATPS)
            resultfile=[OCMATPS.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATPS.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATPS.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATPS'' changed.\n')
        OCMATPS=MODELINFO.OCMATPS;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATPSORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case 'odesolution'
        global OCMATCONT OCMATODEBVP OCBVP
        % store value of the global variables OCMATCONT and OCMATODEBVP to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATODEBVP)
            OCMATODEBVPORIGINAL=OCMATODEBVP;
        else
            OCMATODEBVPORIGINAL=[];
        end
        if ~isempty(OCBVP)
            OCBVPORIGINAL=OCBVP;
        else
            OCBVPORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATODEBVP)
            resultfile=[OCMATODEBVP.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATODEBVP.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATODEBVP.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATODEBVP'' changed.\n')
        OCMATODEBVP=MODELINFO.OCMATODEBVP;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATODEBVPORIGINAL;
        varargout{3}=OCBVPORIGINAL;
        
%%%%%%%%%%
%% solutions calculated by GRADIENT method
%%%%%%%%%

    case {'extremalgrad4ft','indifferencegradsolution','extremalgrad2ep'}
        global OCGRADCONT OCGRADSOL
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCGRADCONT)
            OCGRADCONTORIGINAL=OCGRADCONT;
        else
            OCGRADCONTORIGINAL=[];
        end
        if ~isempty(OCGRADSOL)
            OCGRADSOLORIGINAL=OCGRADSOL;
        else
            OCGRADSOLORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCGRADSOL)
            resultfile=[OCGRADSOL.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCGRADSOL.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has already been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCGRADCONT'' changed.\n')
        OCGRADCONT=MODELINFO.OCGRADCONT;
        ocmatmsg('Global variable ''OCGRADSOL'' changed.\n')
        OCGRADSOL=MODELINFO.OCGRADSOL;
        varargout{1}=OCGRADCONTORIGINAL;
        varargout{2}=OCGRADSOLORIGINAL;
        
        %%%

    case {'extremaldae4ft'}
        global OCMATCONT OCMATFTE
        % store value of the global variables OCMATCONT and OCMATAE to
        % temporary variables, in an intermediate step the global variables
        % are changed to values stored in the continuation files and reset
        % afterwards
        if ~isempty(OCMATCONT)
            OCMATCONTORIGINAL=OCMATCONT;
        else
            OCMATCONTORIGINAL=[];
        end
        if ~isempty(OCMATFTE)
            OCMATFTEORIGINAL=OCMATFTE;
        else
            OCMATFTEORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(OCMATFTE)
            resultfile=[OCMATFTE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[OCMATFTE.basicglobalvarfilename '4' fieldname '.mat'];
        else
            resultfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ocObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
        end

        % test if files exist
        if ~exist(resultfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',resultfile)
            return
        end
        if ~exist(globalvarfile,'file')
            ocmatmsg('Continuation file ''%s'' does not exist\n',globalvarfile)
            return
        end

        % test if data file is older than 12 hours in that case ask if user
        % want to proceed or interrupt
        info=dir(resultfile);
        if abs(info.datenum-datenum(now))>0.5
            ocmatmsg('Data file has been generated at %s.\n',datestr(info.datenum))
            answer=input('Interrupt storing?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Return without storing.')
                return
            end
        end
        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        % change actual global variables by the values stored in the
        % continuation file as MODELINFO.OCMATCONT and MODELINFO.OCMATAE.
        % If the storing process is interrupted before the global variables
        % are reset to its original value the user is informed that s/he
        % has to newly initialize an actual continuation process.
        ocmatmsg('Global variable ''OCMATCONT'' changed.\n')
        OCMATCONT=MODELINFO.OCMATCONT;
        OCMATCONT.meshadaptationflag=1;
        ocmatmsg('Global variable ''OCMATFTE'' changed.\n')
        OCMATFTE=MODELINFO.OCMATFTE;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=OCMATFTEORIGINAL;


end

function resetglobal(fieldname,varargin)
switch fieldname
    case {'extremal2ep','extremalt2ep','extremalp2ep','extremalc2ep','extremalp2epuser','extremale2ep','extremal2per','extremal2lc','extremal2emf','extremalt2emf','optisocline','extremal2inf','extremalt2inf','extremalp2inf','extremal2ftae','saddlepath2ep'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATAEORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATAEORIGINAL)
            ocmatmsg('Global variable ''OCMATAE'' reset.\n')
            OCMATAE=OCMATAEORIGINAL;
        else
            clear global OCMATAE
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end
    case {'extremal4ft','extremalt4ft','extremalp4ft','extremaltF4ft','extremalpF4ft','extremalF4ft','extremalt4inft'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATFTEORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATFTEORIGINAL)
            ocmatmsg('Global variable ''OCMATFTE'' reset.\n')
            OCMATFTE=OCMATFTEORIGINAL;
        else
            clear global OCMATFTE
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'indifferencesolution','indifferencesolution4emf','indifferencesolution4per','indifferencedistribution','indifferencesolution4ft','indifferencesolutionp','indifferencesolution4ae2ftae','indifferencesolution4ep_ft','indifferencesolution4mm', ...
            'indifferencesolution4lcae','indifferencesolution2lc'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATINDIFORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATINDIFORIGINAL)
            ocmatmsg('Global variable ''OCMATINDIF'' reset.\n')
            OCMATINDIF=OCMATINDIFORIGINAL;
        else
            clear global OCMATINDIF
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'heteroclinic','heterocliniccyc','heteroclinicep2ft','homoclinic','heteroclinicep2lc'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATHETORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATHETORIGINAL)
            ocmatmsg('Global variable ''OCMATHET'' reset.\n')
            OCMATHET=OCMATHETORIGINAL;
        else
            clear global OCMATHET
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'extremalmp2ep'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATAEORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATAEORIGINAL)
            ocmatmsg('Global variable ''OCMATAE'' reset.\n')
            OCMATAE=OCMATAEORIGINAL;
        else
            clear global OCMATAE
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'limitextremal','limitextremal','limitextremal2lc','limitextremal4ft'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATLSCORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATLSCORIGINAL)
            ocmatmsg('Global variable ''OCMATLSC'' reset.\n')
            OCMATLSC=OCMATLSCORIGINAL;
        else
            clear global OCMATLSC
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'odesolution'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATODEBVPORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATODEBVPORIGINAL)
            ocmatmsg('Global variable ''OCMATODEBVP'' reset.\n')
            OCMATODEBVP=OCMATODEBVPORIGINAL;
        else
            clear global OCMATODEBVP
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'limitcycle','limitcycleuser','limitcycleuserI','limitlimitcycle'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATLCORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATLCORIGINAL)
            ocmatmsg('Global variable ''OCMATLC'' reset.\n')
            OCMATLC=OCMATLCORIGINAL;
        else
            clear global OCMATLC
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end

    case {'modelperiodicsol','periodicsol_T'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATPSORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATPSORIGINAL)
            ocmatmsg('Global variable ''OCMATPS'' reset.\n')
            OCMATPS=OCMATPSORIGINAL;
        else
            clear global OCMATPS
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end
        
%%%%%%%%%%
%% solutions calculated by GRADIENT method
%%%%%%%%%

    case {'extremalgrad4ft'}
        OCGRADCONTORIGINAL=varargin{1};
        OCGRADSOLORIGINAL=varargin{2};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCGRADCONTORIGINAL)
            ocmatmsg('Global variable ''OCGRADCONT'' reset.\n')
            OCGRADCONT=OCGRADCONTORIGINAL;
        else
            clear global OCGRADCONT
        end
        if ~isempty(OCGRADSOLORIGINAL)
            ocmatmsg('Global variable ''OCGRADSOL'' reset.\n')
            OCGRADSOL=OCGRADSOLORIGINAL;
        else
            clear global OCGRADSOL
        end
        
%%%%%%%%%%
%% solutions calculated by DAE method
%%%%%%%%%

    case {'extremaldae4ft'}
        OCMATCONTORIGINAL=varargin{1};
        OCMATFTEORIGINAL=varargin{2};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(OCMATFTEORIGINAL)
            ocmatmsg('Global variable ''OCMATFTE'' reset.\n')
            OCMATFTE=OCMATFTEORIGINAL;
        else
            clear global OCMATFTE
        end
end

function [solt sole soln sol]=findspecificsolution(sout)
% write violation information into violationinfo field of soln
solt=[]; % solution at target point
sole=[]; % last solution
sol=[];
soln=[];
if isempty(sout)
    return
end
counter=numel(sout);
sole=sout(counter).data.sol;
while counter>2
    counter=counter-1;
    if strcmp(sout(counter).label,'NAS')
        sol=sout(counter).data.sol;
        soln=sout(counter).data.soln;
        sout(counter).data=rmfield(sout(counter).data,{'sol','soln'});
        soln.violationinfo=sout(counter).data.infoS;
    end
    if strcmp(sout(counter).label,'HTV')
        solt=sout(counter).data.sol;
    end
end
