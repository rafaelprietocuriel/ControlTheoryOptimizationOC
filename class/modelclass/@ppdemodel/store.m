function varargout=store(ppdeObj,varargin)
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
%        ppdeasymptotic: asymptotic solutions
%        cell array of 'dynprimitives' or 'ppdeasymptotics'
%
% The default names of the Result field are
%       Result.Equilibrium: for a dynprimitive being an equilibrium
%       Result.LimitCycle: for a dynprimitive being a limit cycle
%       Result.ExtremalSolution: for an ppdeasymptotic
%       Result.Continuation: for the result of a continuation process
% If no output argument is specified it is written into the object OCOBJ.
%
% OCOBJN=STORE(...) the result is returned to the new stdocmodel instance
% OCOBJN 


if isempty(ppdeObj) || nargin==1
    return
end

fieldname=[];
fieldvalue=[];
if ~ischar(varargin{1}) && nargin==2
    if ispdeprimitive(varargin{1}) ||  ispdeasymptotic(varargin{1}) || ispdetrajectory(varargin{1})
        fieldname=class(varargin{1});
    elseif iscell(varargin{1})
        for ii=1:numel(varargin{1})
            store(ppdeObj,varargin{1}{ii},varargin{2:end});
        end

        if ~nargout
            assignin('caller',inputname(1),ppdeObj);
        else
            varargout{1}=ppdeObj;
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
[ocStruct,ocResultFieldName]=generateelement(ppdeObj,fieldname,fieldvalue);
if isempty(ocResultFieldName) || isempty(ocStruct)
    return
end
if isfield(ppdeObj.Result,ocResultFieldName)
    ppdeObj.Result.(ocResultFieldName){end+1}=ocStruct;
else
    ppdeObj.Result(1).(ocResultFieldName){1}=ocStruct;
end

if ~nargout
    assignin('caller',inputname(1),ppdeObj);
else
    varargout{1}=ppdeObj;
end

%%
function [ocResultStruct,ocResultFieldName]=generateelement(ppdeObj,fieldname,fieldvalue)
ocResultStruct=[];
ocResultFieldName='';
switch fieldname
    case 'pdeprimitive'
        if isequilibrium(fieldvalue)
            ocResultStruct=fieldvalue;
            ocResultFieldName='Equilibrium';
        else
            ocResultStruct=fieldvalue;
            ocResultFieldName='PeriodicSolution';
        end
    case 'pdeasymptotic'
        ocResultStruct=fieldvalue;
        ocResultFieldName='ExtremalSolution';
        
    case {'extremal2ee'}
        global OCMATCONT OCMATAE OCBVP 

        [globalvarfile,resultfile,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL]=loadglobal(ppdeObj,fieldname);

        % load data of information files into workspace
        load(globalvarfile)
        load(resultfile)

        mainfunch=feval(str2func(fieldname));
        %mainfunch=feval(str2func('ppdeextremal2epdeS'));

        % determine if non admissible solution occurs
        [soltarget,sollast,soln,sol]=findspecificsolution(sout);
        counter=0;
        for ii=1:numel(bvpout)
            try
                counter=counter+1;
                DataAdaptation(bvpout(ii).tmesh);
                ocExStruct(counter)=mainfunch{22}(bvpout(ii).tmesh,bvpout(ii).coeff,[]);
            catch
                ocmatmsg('Problems generating solution structure from continuation step %d\n',ii)
            end
        end
        if any(strcmp(fieldname,{'extremal2ee'}))
            ocEP=pdeprimitive(OCMATAE.saddlepoint,OCMATCONT.HE.arcarg(end),OCMATAE.femdata);
            ocEP.linearization=OCMATAE.linearization;
        end
       if ~isempty(soln)
            ocResultStruct.NonadmissibleSolution=pdeasymptotic(pdetrajectory(soln),ocEP,soln.solverinfo.pathtype);
            ocResultStruct.ExtremalSolution=pdeasymptotic(pdetrajectory(sol),ocEP,sol.solverinfo.pathtype);
        else
            ocResultStruct.NonadmissibleSolution=pdeasymptotic([]);
        end
        if ~isempty(soltarget)
            ocResultStruct.ExtremalSolution=pdeasymptotic(pdetrajectory(soltarget),ocEP,soltarget.solverinfo.pathtype);
        else
            ocResultStruct.ExtremalSolution=pdeasymptotic(pdetrajectory(ocExStruct(end)),ocEP,ocExStruct(end).solverinfo.pathtype);
        end
        ocResultStruct.ContinuationSolution=ocExStruct;
        ocResultStruct.ContinuationClassification=fieldname;
        ocResultStruct.ContinuationInformation=sout;
        ocResultStruct.LimitSet=ocEP;
        ocResultFieldName='Continuation';

        resetglobal(fieldname,OCMATCONTORIGINAL,OCMATAEORIGINAL,OCBVPORIGINAL);
    otherwise
        ocResultStruct=fieldvalue;
        ocResultFieldName=fieldname;
end

function [globalvarfile,resultfile,varargout]=loadglobal(ppdeObj,fieldname)
switch fieldname
    case {'extremal2ee'}
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
            resultfile=fullfile(getocmatpath,getocmatfolder(ppdeObj,'userdata'),['SaveIntermediateResults4' fieldname '.mat']);
        end
        if isempty(OCMATAE) || ~exist(globalvarfile,'file')
            globalvarfile=fullfile(getocmatpath,getocmatfolder(ppdeObj,'userdata'),['SaveIntermediateResultsGlobalVariable4' fieldname '.mat']);
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
end

function resetglobal(fieldname,varargin)
switch fieldname
    case {'extremal2ee'}
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
