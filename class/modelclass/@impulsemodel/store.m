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
    if isdynprimitive(varargin{1}) ||  isocasymptotic(varargin{1}) || isoctrajectory(varargin{1}) || isoccurve(varargin{1}) ||  isocmultipath(varargin{1})
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

    case {'extremal4fimp','extremalt4fimp','extremalp4fimp','extremalpopt4fimp'}
        global OCMATCONT IOCMATFTE OCBVP
        [globalvarfile,resultfile,OCMATCONTORIGINAL,IOCMATFTEORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

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
                slMfStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                slMfStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                slMfStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                endpointStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(end)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(end)).eqcoord,end);
                endpointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(end);
                endpointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
                if IOCMATFTE.objectivevaluecalc
                    slMfStruct.userinfo.objectivevalue(counter)=ocExStruct(ii).y(IOCMATFTE.objectivevaluecoord,end);
                end
            end
        end
        %OCMATCONT.meshadaptationflag=meshadaptationflag;
        slMfStruct.arcposition=[];
        endpointStruct.arcposition=[];
        if ~isempty(soln)
            if strcmp(fieldname,'extremalp4fimp')
                ocObjtmp=changeparametervalue(ocObj,soln.modelparameter);
            elseif strcmp(fieldname,'extremalpopt4fimp')
                ocObjtmp=changeparametervalue(ocObj,IOCMATFTE.optparindex,soln.solverinfo.parameters(IOCMATFTE.optparametercoord));
                if strncmp(IOCMATFTE.targettype,'p',1)
                    ocObjtmp=changeparametervalue(ocObj,IOCMATFTE.contindex,soln.solverinfo.parameters(IOCMATFTE.contcoord));
                end
            else
                ocObjtmp=ocObj;
            end
            soln.jumparg(sol.jumparg<0)=ceil(soln.jumparg(sol.jumparg<0));
            ocResultStruct.NonadmissibleSolution=hybridoctrajectory(sol,ocObjtmp);
            if strcmp(fieldname,'extremalp4fimp')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            elseif strcmp(fieldname,'extremalpopt4fimp')
                ocObjtmp=changeparametervalue(ocObj,IOCMATFTE.optparindex,sol.solverinfo.parameters(IOCMATFTE.optparametercoord));
                if strncmp(IOCMATFTE.targettype,'p',1)
                    ocObjtmp=changeparametervalue(ocObj,IOCMATFTE.contindex,sol.solverinfo.parameters(IOCMATFTE.contcoord));
                end
            else
                ocObjtmp=ocObj;
            end
            sol.jumparg(sol.jumparg<0)=ceil(sol.jumparg(sol.jumparg<0));
            ocResultStruct.ExtremalSolution=hybridoctrajectory(sol,ocObjtmp);
        else
            if ~isempty(soltarget)
                sol=soltarget;
            else
                sol=ocExStruct(end);
            end
            sol.jumparg(sol.jumparg<0)=ceil(sol.jumparg(sol.jumparg<0));
            if strcmp(fieldname,'extremalp4fimp')
                ocObjtmp=changeparametervalue(ocObj,sol.modelparameter);
            elseif strcmp(fieldname,'extremalpopt4fimp')
                ocObjtmp=changeparametervalue(ocObj,IOCMATFTE.optparindex,sol.solverinfo.parameters(IOCMATFTE.optparametercoord));
                if strncmp(IOCMATFTE.targettype,'p',1)
                    ocObjtmp=changeparametervalue(ocObj,IOCMATFTE.contindex,sol.solverinfo.parameters(IOCMATFTE.contcoord));
                end
            else
                ocObjtmp=ocObj;
            end
            ocResultStruct.ExtremalSolution=hybridoctrajectory(sol,ocObjtmp);
            ocResultStruct.NonadmissibleSolution=hybridoctrajectory([]);
        end
        ocResultStruct.SliceManifold=occurve(slMfStruct);
        ocResultStruct.EndpointManifold=occurve(endpointStruct);
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        if strcmp(fieldname,'extremalp4fimp')
            ocResultStruct.ContinuationParameter=IOCMATFTE.varyparameterindex;
            ocResultStruct.LastModel=ocObjtmp;
        elseif strcmp(fieldname,'extremalpopt4fimp')
            ocResultStruct.ContinuationParameter=IOCMATFTE.optparametercoord;
            if strncmp(IOCMATFTE.targettype,'p',1)
                ocResultStruct.ContinuationParameter=[ocResultStruct.ContinuationParameter IOCMATFTE.contindex];
            end
            ocResultStruct.LastModel=ocObjtmp;
        end
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,IOCMATFTEORIGINAL,OCBVPORIGINAL);

    case {'indifferencesolution4fimp'}
        global OCMATCONT IOCMATINDIF OCBVP

        [globalvarfile,resultfile,OCMATCONTORIGINAL,IOCMATINDIFORIGINAL,OCBVPORIGINAL]=loadglobal(ocObj,fieldname);

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
                indifferencepointStruct.y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,counter)=ocExStruct(ii).y(OCMATCONT.DOMAINDDATA(OCMATCONT.HE.arcindex(1)).eqcoord,1);
                indifferencepointStruct.arcarg(counter)=OCMATCONT.HE.arcarg(1);
                indifferencepointStruct.userinfo.continuationparameter(counter)=bvpout(ii).coeff(end);
            end
        end
        indifferencepointStruct.arcposition=[];

        if ~isempty(soln)
            if any(strcmp(fieldname,{'indifferencesolution4fimp'}))
                ocResultStruct.NonadmissibleSolution=sol2hybridocmultipath(soln,IOCMATINDIF);
                ocResultStruct.ExtremalSolution=sol2hybridocmultipath(sol,IOCMATINDIF);
            end
        else
            if ~isempty(soltarget)
                if any(strcmp(fieldname,{'indifferencesolution4fimp'}))
                    ocResultStruct.ExtremalSolution=sol2hybridocmultipath(soltarget,IOCMATINDIF);
                end
            else

                if any(strcmp(fieldname,{'indifferencesolution4fimp'}))
                    ocResultStruct.ExtremalSolution=sol2hybridocmultipath(ocExStruct(end),IOCMATINDIF);
                end
            end
            ocResultStruct.NonadmissibleSolution=hybridoctrajectory([]);
        end
        if any(strcmp(fieldname,{'indifferencesolution4fimp'}))
            ocResultStruct.IndifferencePointCurve=occurve(indifferencepointStruct);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,IOCMATINDIFORIGINAL,OCBVPORIGINAL);
    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(ocObj,fieldname)
switch fieldname

    case {'extremal4fimp','extremalt4fimp','extremalp4fimp','extremalpopt4fimp'}
        global OCMATCONT OCBVP IOCMATFTE
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
        if ~isempty(IOCMATFTE)
            IOCMATFTEORIGINAL=IOCMATFTE;
        else
            IOCMATFTEORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(IOCMATFTE)
            resultfile=[IOCMATFTE.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[IOCMATFTE.basicglobalvarfilename '4' fieldname '.mat'];
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
        ocmatmsg('Global variable ''IOCMATFTE'' changed.\n')
        IOCMATFTE=MODELINFO.IOCMATFTE;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=IOCMATFTEORIGINAL;
        varargout{3}=OCBVPORIGINAL;

    case {'indifferencesolution4fimp'}
        global OCMATCONT OCBVP IOCMATINDIF
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
        if ~isempty(IOCMATINDIF)
            IOCMATINDIFORIGINAL=IOCMATINDIF;
        else
            IOCMATINDIFORIGINAL=[];
        end

        % the file names containing continuation information are build by
        % the model specific basic names and the continuation
        % classification, e.g. 4extremal2ep
        if ~isempty(IOCMATINDIF)
            resultfile=[IOCMATINDIF.basicresultfilename '4' fieldname '.mat'];
            globalvarfile=[IOCMATINDIF.basicglobalvarfilename '4' fieldname '.mat'];
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
        ocmatmsg('Global variable ''IOCMATINDIF'' changed.\n')
        IOCMATINDIF=MODELINFO.IOCMATINDIF;
        ocmatmsg('Global variable ''OCBVP'' changed.\n')
        OCBVP=MODELINFO.OCBVP;
        varargout{1}=OCMATCONTORIGINAL;
        varargout{2}=IOCMATINDIFORIGINAL;
        varargout{3}=OCBVPORIGINAL;

end

function resetglobal(fieldname,varargin)
switch fieldname
    case {'extremal4fimp','extremalt4fimp','extremalp4fimp','extremalpopt4fimp'}
        OCMATCONTORIGINAL=varargin{1};
        IOCMATFTEORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(IOCMATFTEORIGINAL)
            ocmatmsg('Global variable ''IOCMATFTE'' reset.\n')
            IOCMATFTE=IOCMATFTEORIGINAL;
        else
            clear global IOCMATFTE
        end
        if ~isempty(OCBVPORIGINAL)
            ocmatmsg('Global variable ''OCBVP'' reset.\n')
            OCBVP=OCBVPORIGINAL;
        else
            clear global OCBVP
        end
    case {'indifferencesolution4fimp'}
        OCMATCONTORIGINAL=varargin{1};
        IOCMATINDIFORIGINAL=varargin{2};
        OCBVPORIGINAL=varargin{3};
        % reset global variables to its original value or delete it if no
        % global variables existed
        if ~isempty(OCMATCONTORIGINAL)
            ocmatmsg('Global variable ''OCMATCONT'' reset.\n')
            OCMATCONT=OCMATCONTORIGINAL;
        else
            clear global OCMATCONT
        end
        if ~isempty(IOCMATINDIFORIGINAL)
            ocmatmsg('Global variable ''IOCMATINDIF'' reset.\n')
            IOCMATINDIF=IOCMATINDIFORIGINAL;
        else
            clear global IOCMATINDIF
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
