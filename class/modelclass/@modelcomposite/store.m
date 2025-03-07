function varargout=store(mmObj,varargin)
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


if isempty(mmObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1}) || isoctrajectory(varargin{1}) || isoccurve(varargin{1}) ||  isocmultipath(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(mmObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),mmObj);
        else
            varargout{1}=mmObj;
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
[ocStruct,ocResultFieldName]=generateelement(mmObj,fieldname,fieldvalue);
if isempty(ocResultFieldName) || isempty(ocStruct)
    return
end
if isfield(mmObj.Result,ocResultFieldName)
    mmObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    mmObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),mmObj);
else
    varargout{1}=mmObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(mmObj,fieldname,fieldvalue)
ocResultStruct=[];
ocResultFieldName='';
switch fieldname


    case {'extremal4ft_mm'}
        global OCMATCONT OCMATFTE OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATFTEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);

        counter=0;
        %meshadaptationflag=OCMATCONT.meshadaptationflag;
        OCMATCONT.meshadaptationflag=1;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(OCMATCONT.DOMAINDDATA{1}(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA{1}(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                endpointStruct.y(OCMATCONT.DOMAINDDATA{1}(OCMATCONT.HE.arcindex(end)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA{1}(OCMATCONT.HE.arcindex(end)).eqcoord,end);
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
            ocResultStruct.ExtremalSolution=sol2mmultipath(sol,mmObj);
            ocResultStruct.NonadmissibleSolution=sol2mmultipath([]);
        elseif ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=sol2mmultipath(soln,mmObj);
            ocResultStruct.ExtremalSolution=sol2mmultipath(sol,mmObj);
        else
            sol=ocExStruct(end);
            ocResultStruct.ExtremalSolution=sol2mmultipath(sol,mmObj);
            ocResultStruct.NonadmissibleSolution=sol2mmultipath([]);
        end
        ocResultStruct.SliceManifold=occurve(slMfStruct);
        ocResultStruct.EndpointManifold=occurve(endpointStruct);
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATFTEORIGINAL,OCBVPORIGINAL);


    case {'indifferencesolution4mm'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

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

        %         ocEP=dynprimitive(OCMATINDIF.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATINDIF.linearization);
%         if ~isempty(soln)
%             if strcmp(fieldname,'indifferencesolution4mm')
%                 ocResultStruct.NonadmissibleSolution=sol2occomposite(soln,mmObj);
%                 ocResultStruct.ExtremalSolution=sol2occomposite(sol,mmObj);
%             end
%         else
%             if ~isempty(soltarget)
%                 if strcmp(fieldname,'indifferencesolution4mm')
%                     ocResultStruct.ExtremalSolution=sol2occomposite(soltarget,mmObj);
%                 end
%             else
%                 sol=ocExStruct(end);
%                 if strcmp(fieldname,'indifferencesolution4mm')
%                     ocResultStruct.ExtremalSolution=sol2occomposite(sol,mmObj);
%                 end
%             end
%         end
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


    case {'indifferencesolutionp'}
        global OCMATCONT OCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));


        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
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
            mmObjtmp=changeparametervalue(mmObj,OCMATINDIF.parameterindex,soln.solverinfo.parameters(end));
            par=parametervalue(mmObjtmp);
            J=cell(1,OCMATINDIF.indifferenceorder);
            saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
            for ii=1:OCMATINDIF.indifferenceorder
                equcoord=soln.solverinfo.equilibriumcoord{ii};
                hatx.y=soln.solverinfo.parameters(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                J{ii}=OCMATINDIF.canonicalsystemjacobian(0,hatx.y,par,hatx.arcarg);
                if ~isempty(OCMATINDIF.excludecoordinate4ep{ii})
                    J{ii}(:,OCMATINDIF.excludecoordinate4ep{ii})=[];
                    J{ii}(OCMATINDIF.excludecoordinate4ep{ii},:)=[];
                end
                saddlepoint{ii}=hatx.y;
            end
            OCMATINDIF.indifferenceorder=OCMATINDIF.indifferenceorder;
            ocResultStruct.NonadmissibleSolution=sol2mmultipath(soln,OCMATINDIF,saddlepoint,J);

            mmObjtmp=changeparametervalue(mmObj,OCMATINDIF.parameterindex,sol.solverinfo.parameters(end));
            par=parametervalue(mmObjtmp);
            J=cell(1,OCMATINDIF.indifferenceorder);
            saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
            for ii=1:OCMATINDIF.indifferenceorder
                equcoord=sol.solverinfo.equilibriumcoord{ii};
                hatx.y=sol.solverinfo.parameters(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                J{ii}=OCMATINDIF.canonicalsystemjacobian(0,hatx.y,par,hatx.arcarg);
                if ~isempty(OCMATINDIF.excludecoordinate4ep{ii})
                    J{ii}(:,OCMATINDIF.excludecoordinate4ep{ii})=[];
                    J{ii}(OCMATINDIF.excludecoordinate4ep{ii},:)=[];
                end
                saddlepoint{ii}=hatx.y;
            end
            ocResultStruct.ExtremalSolution=sol2mmultipath(sol,OCMATINDIF,saddlepoint,J);
        else
            mmObjtmp=changeparametervalue(mmObj,OCMATINDIF.parameterindex,soltarget.solverinfo.parameters(end));
            par=parametervalue(mmObjtmp);
            J=cell(1,OCMATINDIF.indifferenceorder);
            saddlepoint=cell(1,OCMATINDIF.indifferenceorder);
            for ii=1:OCMATINDIF.indifferenceorder
                equcoord=soltarget.solverinfo.equilibriumcoord{ii};
                hatx.y=soltarget.solverinfo.parameters(equcoord);
                hatx.x=0;
                hatx.arcarg=OCMATCONT.HE.arcarg(OCMATINDIF.arccoord{ii}(end));
                J{ii}=OCMATINDIF.canonicalsystemjacobian(0,hatx.y,par,hatx.arcarg);
                if ~isempty(OCMATINDIF.excludecoordinate4ep{ii})
                    J{ii}(:,OCMATINDIF.excludecoordinate4ep{ii})=[];
                    J{ii}(OCMATINDIF.excludecoordinate4ep{ii},:)=[];
                end
                saddlepoint{ii}=hatx.y;
            end
            OCMATINDIF.indifferenceorder=OCMATINDIF.indifferenceorder;
            ocResultStruct.ExtremalSolution=sol2mmultipath(soltarget,OCMATINDIF,saddlepoint,J);
            ocResultStruct.NonadmissibleSolution=ocasymptotic([]);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LastModel=mmObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATINDIFORIGINAL,OCBVPORIGINAL);

    case {'heteroclinic','heterocliniccyc','heteroclinicep2ft','homoclinic'}
        global OCMATCONT OCMATHET OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATHETORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        %OCMATHET.hetorder=2;
        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
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
            ocResultStruct.NonadmissibleSolution=sol2mmultipath(soln,OCMATHET,OCMATHET.saddlepoint,OCMATHET.linearization);
            ocResultStruct.ExtremalSolution=sol2mmultipath(sol,OCMATHET,OCMATHET.saddlepoint,OCMATHET.linearization);
        else
            parnew=soltarget.modelparameter;
            parnew(OCMATHET.parameterindex)=soltarget.solverinfo.parameters(OCMATHET.parametervaluecoord);
            mmObjtmp=changeparametervalue(mmObj,parnew);
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
                    J{ii}=jacobian(mmObjtmp,hatx);
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
        ocResultStruct.LastModel=mmObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATHETORIGINAL,OCBVPORIGINAL);


    case {'extremalmp2ep'}
        global OCMATCONT OCMATAE OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        OCMATAE.hetorder=2;
        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
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
            mmObjtmp=changeparametervalue(mmObj,parnew);
            J=cell(1,OCMATAE.multorder);
            saddlepoint=cell(1,OCMATAE.multorder);
            for ii=1:OCMATAE.multorder
                calcequilibrium=1;
                equcoord=soltarget.solverinfo.equilibriumcoord;
                if calcequilibrium
                    hatx.y=soltarget.solverinfo.parameters(equcoord);
                    hatx.x=0;
                    hatx.arcarg=OCMATCONT.HE.arcarg(OCMATAE.arccoord{ii}(end));
                    J{ii}=jacobian(mmObjtmp,hatx);
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
        ocResultStruct.LastModel=mmObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'limitcycle','limitcycleuser','limitcycleuserI'}
        global OCMATCONT OCMATLC OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
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
            if any(strcmp(fieldname,{'limitcycle','limitcycleuser'}))
                mmObjtmp=changeparametervalue(mmObj,soln.octrajectory.modelparameter);
            else
                mmObjtmp=mmObj;
            end
            ocResultStruct.NonadmissibleLimitCycle=dynprimitive(soln,mmObjtmp);
            ocResultStruct.LimitCycle=dynprimitive(sol,mmObjtmp);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.LimitCycleMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                if any(strcmp(fieldname,{'limitcycle','limitcycleuser'}))
                    mmObjtmp=changeparametervalue(mmObj,soltarget.octrajectory.modelparameter);
                else
                    mmObjtmp=mmObj;
                end
                ocResultStruct.LimitCycle=dynprimitive(soltarget,mmObjtmp);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.LimitCycleMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                if any(strcmp(fieldname,{'limitcycle','limitcycleuser'}))
                    mmObjtmp=changeparametervalue(mmObj,ocExStruct(end).octrajectory.modelparameter);
                else
                    mmObjtmp=mmObj;
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
        ocResultStruct.LastModel=mmObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'modelperiodicsol'}
        global OCMATCONT OCMATPS OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if ~isempty(soln)
            mmObjtmp=changeparametervalue(mmObj,OCMATPS.varyparameterindex,soln.octrajectory.solverinfo.continuationparameter);
            ocResultStruct.NonadmissiblePeriodicSolution=dynprimitive(soln,mmObjtmp);
            ocResultStruct.PeriodicSolution=dynprimitive(sol,mmObjtmp);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                mmObjtmp=changeparametervalue(mmObj,OCMATPS.varyparameterindex,soltarget.octrajectory.solverinfo.continuationparameter);
                ocResultStruct.PeriodicSolution=dynprimitive(soltarget,mmObjtmp);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
                end
            else
                %                 DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                mmObjtmp=changeparametervalue(mmObj,OCMATPS.varyparameterindex,ocExStruct(end).octrajectory.solverinfo.continuationparameter);
                %                 ocResultStruct.PeriodicSolution=dynprimitive(soltarget,mmObjtmp);
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
        ocResultStruct.LastModel=mmObjtmp;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'periodicsol_T'}
        global OCMATCONT OCMATPS OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
            end
        end
        if ~isempty(soln)
            ocResultStruct.NonadmissiblePeriodicSolution=dynprimitive(soln,mmObj);
            ocResultStruct.PeriodicSolution=dynprimitive(sol,mmObj);
            trajectory=ocExStruct(end).octrajectory;
            matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocResultStruct.PeriodicSolution=dynprimitive(soltarget,mmObj);
                trajectory=soltarget.octrajectory;
                matlabsol=transform2nativematlab(trajectory.solverinfo.tmesh,trajectory.solverinfo.coeff,trajectory.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.PeriodicSolutionMatlabFormat=matlabsol;
                end
            else
                %                 DataAdaptation(ocExStruct(end).octrajectory.solverinfo.tmesh);
                %                 ocResultStruct.PeriodicSolution=dynprimitive(soltarget,mmObjtmp);
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

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,bvpout(ii).tangent);
                slMfStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        slMfStruct.arcposition=[];
        if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=ocasymptotic(octrajectory(soln,mmObj),OCMATAE.limitset);
            ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(sol,mmObj),OCMATAE.limitset);
            matlabsol=transform2nativematlab(ocExStruct(end).x,ocExStruct(end).solverinfo.coeff,ocExStruct(end).modelparameter);
            if ~isempty(matlabsol)
                ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
            end
        else
            if ~isempty(soltarget)
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(soltarget,mmObj),OCMATAE.limitset);
                matlabsol=transform2nativematlab(soltarget.solverinfo.tmesh,soltarget.solverinfo.coeff,soltarget.modelparameter);
                if ~isempty(matlabsol)
                    ocResultStruct.ExtremalSolutionMatlabFormat=matlabsol;
                end
            else
                DataAdaptation(ocExStruct(end).solverinfo.tmesh);
                ocResultStruct.ExtremalSolution=ocasymptotic(octrajectory(ocExStruct(end),mmObj),OCMATAE.limitset);
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
        ocResultStruct.LimitSet=OCMATAE.limitset;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);

    case {'odesolution'}
        global OCMATCONT OCMATODEBVP OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
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

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(mmObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));

        % determine if non admissible solution occurs
        [soltarget sollast soln sol]=findspecificsolution(sout);
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
    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(mmObj,fieldname)
switch fieldname
    case {'extremal2ep','extremalt2ep','extremalp2ep','extremalc2ep','extremalp2epuser','extremale2ep','extremal2per','extremal2lc','extremal2emf','extremalt2emf','optisocline','extremal2inf','extremalt2inf','extremalp2inf','extremal2ftae'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
        end
        if isempty(OCMATAE) || ~exist(globalvarfile,'file')
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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

    case {'extremal4ft_mm'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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
            funch=feval(modelspecificfunc(mmObj,'4SaddlePathContinuation'));
        else
            funch=feval(modelspecificfunc(mmObj,'4FiniteHorizonPathContinuation'));
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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

    case {'indifferencesolution4mm'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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


    case {'heteroclinic','heterocliniccyc','heteroclinicep2ft','homoclinic'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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

    case {'limitcycle','limitcycleuser','limitcycleuserI'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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
            resultfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
            globalvarfile=fullfile(getocmatpath,getocmatfolder(mmObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
            ocmatmsg('Data file was already generated at %s.\n',datestr(info.datenum))
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
end

function resetglobal(fieldname,varargin)
switch fieldname
    case {'extremal2ep','extremalt2ep','extremalp2ep','extremalc2ep','extremalp2epuser','extremale2ep','extremal2per','extremal2lc','extremal2emf','extremalt2emf','optisocline','extremal2inf','extremalt2inf','extremalp2inf','extremal2ftae'}
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
    case {'extremal4ft_mm'}
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

    case {'indifferencesolution4mm'}
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

    case {'heteroclinic','heterocliniccyc','heteroclinicep2ft','homoclinic'}
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

    case {'limitcycle','limitcycleuser','limitcycleuserI'}
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
