function sol=initocmat_MM_FTE_IS(ocObj,compTrj,continuationtype,continuationindex,continuationtarget,varargin)
%
% initocmat_MM_FTE_IS initialization for indifference solution continuation
% of multimodel solution
%
% SOL=initocmat_MM_FTE_IS(MMOBJ,OCTRJM,CONTINUATIONTYPE,CONTINUATIONINDEX,CONTINUATIONTARGET) the
% continuation process is initialized (started) at a solution OCTRJM (cell
% of mmultipath solutions
%
% We assume that the different solutions are of the same structure of parts
% (stages) and the according models are the same for each corresponding
% part
% EXAMPLE: three parts, indifference order two: if the first solution is a
% solution of the models m_1/par_1, m_2/par_2, m_3/par_3, then the second
% solution is of the same sturcture, same models with same parameter values
%
% We will also allow one part/stage solutions of any indifference order:
% this means that for every solution we allow a different model and or
% different parameter values
%
% FREEPARAMETERVECTOR: since the parameter values or/an its location can be
% different for every part in every indifference solution, we have to
% specify the indices of free parameters for each parametervector.
% Example: two parts, indifference order two: alpha,beta appear in all
% parts at the same location, say 1,2 then the freeparametervector is given
% by a cellarray of length four: {[1 2]},{[1 2]},{[1 2]},{[1 2]} the target
% value is of the same structure
%
%two parts, indifference order two: alpha,beta appear in different parts at
%a different location´, say part1: [1 2], part2: [2 1], then the
%freeparametervector is given by a cellaray of length four: {[1 2]},{[2
%1]},{[1 2]},{[2 1]}
%
%two parts, indifference order two: in the first part we want to continue
%along alpha, beta, [1 2] in the second part along sigma and gamma [3 4]
%then the freeparametervector is given by a cellarray of length four:  
% {[1 2]},{[3 4]},{[1 2]},{[3 4]} the target 
% value is of the same structure
% For a more detailed discussion see INITOCMAT_AE_EP
%
clear global OCMATCONT OCMATINDIF
global OCMATCONT OCMATINDIF
sol=[];

fixendstate=[];
fixinitialstate=[];
freeparametervector=[];
freestatevector=[];
freetimevector=[];
freetargetvalue=[];
connectiontimeindex=[];
objectivevaluecalc=[];
objectivevaluecoordinate=[];
exogenousfunction=[];
userbc=[];
if isempty(ocObj)
    ocmatmsg('oc model is empty.')
    return
end

% contiuationtype: time, state, parameter
% continuationindex: (time) number of switching/endtime, (state)
% coordinate, (parameter) index or name of parameter
% continuationtarget: value
fixinitialstateidx=find(strcmpi(varargin,'fixinitialstate'));
fixendstateidx=find(strcmpi(varargin,'fixendstate'));
freeparametervectoridx=find(strcmpi(varargin,'freeparametervector'));
freestatevectoridx=find(strcmpi(varargin,'freestatevector'));
freetimevectoridx=find(strcmpi(varargin,'freetimevector'));
freetargetvalueidx=find(strcmpi(varargin,'freetargetvalue'));
objectivevaluecalcidx=find(strcmpi(varargin,'objectivevaluecalc'));
objectivevaluecoordinateidx=find(strcmpi(varargin,'objectivevaluecoordinate'));
explicitconnectimeidx=find(strcmpi(varargin,'explicitconnectime'));
exogenousfunctionidx=find(strcmpi(varargin,'exogenousfunction'));
connectiontimeindexidx=find(strcmpi(varargin,'connectiontimeindex'));
exogenousinitialstatesidx=find(strcmpi(varargin,'exogenousinitialstates'));
userbcidx=find(strcmpi(varargin,'userbc'));
variationalcalculationidx=find(strcmpi(varargin,'variationalcalculation'));
vjumpargumentidx=find(strcmpi(varargin,'variationaljumpargument'));
vfreetimeidx=find(strcmpi(varargin,'variationalfreetime'));
includevariationalobjectivevalueidx=find(strcmpi(varargin,'includevariationalobjectivevalue'));
if ~isempty(fixinitialstateidx)
    fixinitialstate=varargin{fixinitialstateidx+1};
end
if ~isempty(fixendstateidx)
    fixendstate=varargin{fixendstateidx+1};
end
if ~isempty(freeparametervectoridx)
    freeparametervector=varargin{freeparametervectoridx+1};
end
if ~isempty(freestatevectoridx)
    freestatevector=varargin{freestatevectoridx+1};
end
if ~isempty(freetimevectoridx)
    freetimevector=varargin{freetimevectoridx+1};
end
if ~isempty(freetargetvalueidx)
    freetargetvalue=varargin{freetargetvalueidx+1};
end
if ~isempty(objectivevaluecalcidx)
    objectivevaluecalc=varargin{objectivevaluecalcidx+1};
end
if ~isempty(objectivevaluecoordinateidx)
    objectivevaluecoordinate=varargin{objectivevaluecoordinateidx+1};
end
if ~isempty(explicitconnectimeidx)
    explicitconnectime=varargin{explicitconnectimeidx+1};
end
if ~isempty(exogenousfunctionidx)
    exogenousfunction=varargin{exogenousfunctionidx+1};
end
if ~isempty(connectiontimeindexidx)
    connectiontimeindex=varargin{connectiontimeindexidx+1};
end
if ~isempty(exogenousinitialstatesidx)
    exogenousinitialstates=varargin{exogenousinitialstatesidx+1};
end
if ~isempty(userbcidx)
    userbc=varargin{userbcidx+1};
end
if ~isempty(variationalcalculationidx)
    variationalcalculation=varargin{variationalcalculationidx+1};
end
if ~isempty(vjumpargumentidx)
    vjumpargument=varargin{vjumpargumentidx+1};
end
if ~isempty(vfreetimeidx)
    vfreetime=varargin{vfreetimeidx+1};
end
if ~isempty(includevariationalobjectivevalueidx)
    includevariationalobjectivevalue=varargin{includevariationalobjectivevalueidx+1};
end

if isempty(objectivevaluecalc)
    objectivevaluecalc=0;
end
if iscell(compTrj)
    compTrj=composite(compTrj);
end

if order(ocObj)==order(compTrj)
    orderofmodel=order(ocObj);
    indifferenceorder=orderofmodel;
elseif ismultimodel(ocObj)
    indifferenceorder=order(compTrj);
    orderofmodel=1;
else
    return
end
if isempty(indifferenceorder)
    return
end
OCMATINDIF.exogenousfunction=exogenousfunction;
OCMATINDIF.freetargetvalue=freetargetvalue;

OCMATCONT.modelname=modelname(ocObj);
OCMATINDIF.dhamiltoniandctcoord=[];
counter=0;
totalnumberofstages=0;
connectiontimenum=-inf;
for ii=1:indifferenceorder
    OCMATINDIF.explicitconnectiontime=any(explicitconnectiontime(ocObj(ii)));
    connectiontimenum=max([connectiontimenum,connectiontimenumber(ocObj(ii))]);
    numberofstages=numberofparts(compTrj(ii));
    nummodel(ii)=numberofmodels(ocObj(ii));
    for jj=1:numberofstages
        if orderofmodel==1
            if ii==1
                OCMATCONT.modelfunc{jj}=modelspecificfunc(ocObj(jj),'4IndifferenceSolutionContinuation');
                funch{jj}=OCMATCONT.modelfunc{1}{jj}(); % model specific function handles for saddle path continuation
            end
            modelparameterindex{ii}=1:numberofstages;
        else
            if numberofmodels(ocObj(ii))==numberofstages
                counter=counter+1;
                if ismultistagemodel(ocObj(ii))
                    OCMATCONT.modelfunc{counter}=modelspecificfunc(ocObj(ii),'4FiniteHorizonPathContinuation');
                else
                    OCMATCONT.modelfunc{counter}=modelspecificfunc(ocObj(ii,jj),'4FiniteHorizonPathContinuation');
                end
                funch{counter}=OCMATCONT.modelfunc{counter}(); % model specific function handles for saddle path continuation
            else
                return
            end
        end
    end
    if orderofmodel>1
        totalnumberofstages0=totalnumberofstages+1;
        totalnumberofstages=totalnumberofstages+numberofstages;
        modelparameterindex{ii}=totalnumberofstages0:totalnumberofstages;
    end
end
OCMATINDIF.connectiontimenumber=connectiontimenum;
OCMATINDIF.orderofmodel=orderofmodel;
OCMATINDIF.indifferenceorder=indifferenceorder;
OCMATINDIF.totalnumbermodel=indifferenceorder*orderofmodel;
OCMATINDIF.modelparameterindex=modelparameterindex;
OCMATINDIF.userbc=userbc;
OCMATINDIF.includevariationalobjectivevalue=includevariationalobjectivevalue;

if orderofmodel==1
    forindex=1;
else
    forindex=1:indifferenceorder;
end
counter=0;
for ii=forindex
    numberofstages=numberofparts(compTrj(ii));
    for jj=1:numberofstages
        counter=counter+1;
        OCMATINDIF.canonicalsystem{counter}=funch{counter}{1}{1};
        OCMATINDIF.dhamiltoniandct{counter}=funch{counter}{1}{2};
        OCMATINDIF.dsalvagedct{counter}=funch{counter}{1}{3};

        OCMATINDIF.canonicalsystemjacobian{counter}=funch{counter}{2}{1};
        OCMATINDIF.canonicalsystemparameterjacobian{counter}=funch{counter}{2}{2};
        OCMATINDIF.derivativeconnectiontime{counter}=funch{counter}{2}{4};
        OCMATINDIF.dhamiltoniandctjacobian{counter}=funch{counter}{2}{5};
        OCMATINDIF.dhamiltoniandctparameterjacobian{counter}=funch{counter}{2}{6};
        OCMATINDIF.d2hamiltoniandct2{counter}=funch{counter}{2}{7};

        OCMATINDIF.canonicalsystemhessian{counter}=funch{counter}{3}{1};
        OCMATINDIF.canonicalsystemparameterhessian{counter}=funch{counter}{3}{2};

        % function for the boundary conditions
        OCMATINDIF.bcinitial{counter}=funch{counter}{5}{1};
        OCMATINDIF.bctransversality{counter}=funch{counter}{5}{2};
        OCMATINDIF.bcoptimalhorizon{counter}=funch{counter}{5}{3};
        OCMATINDIF.hitvaluefunc{counter}=funch{counter}{5}{5};
        OCMATINDIF.salvagevalue{counter}=funch{counter}{5}{6};
        OCMATINDIF.bcconnectingparts{counter}=funch{counter}{5}{7};
        OCMATINDIF.bcoptimalconnectingparts{counter}=funch{counter}{5}{8};

        % function for Jacobian
        OCMATINDIF.bcjacobianinitial{counter}=funch{counter}{6}{1};
        OCMATINDIF.bcjacobiantransversality{counter}=funch{counter}{6}{2};

        % function describing the hybrid structure of the problem
        OCMATINDIF.hybridinfo{counter}=funch{counter}{7}{1};
        OCMATINDIF.domain{counter}=funch{counter}{7}{2};
        OCMATINDIF.guard{counter}=funch{counter}{7}{3};
        OCMATINDIF.reset{counter}=funch{counter}{7}{4};
        OCMATINDIF.switchtime{counter}=funch{counter}{7}{5};
        OCMATINDIF.jacobianguard{counter}=funch{counter}{7}{7};
        OCMATINDIF.jacobianreset{counter}=funch{counter}{7}{8};
        OCMATINDIF.domaindiscretization{counter}=funch{counter}{7}{9};

        OCMATINDIF.objectivefunction{counter}=funch{counter}{8}{1};
        OCMATINDIF.objectivefunctionjacobian{counter}=funch{counter}{8}{2};
        OCMATINDIF.objectivefunctionparameterjacobian{counter}=funch{counter}{8}{3};
        OCMATINDIF.objectivefunctionderivativetime{counter}=funch{counter}{8}{4};
        OCMATINDIF.objectivefunctionderivativeconnectiontime{counter}=funch{counter}{8}{5};
        if OCMATINDIF.includevariationalobjectivevalue(counter)
            OCMATINDIF.variationalobjectivefunction{counter}=funch{counter}{8}{5};
            OCMATINDIF.variationalobjectivefunctionjacobian{counter}=funch{counter}{8}{6};
            OCMATINDIF.variationalobjectivefunctionparameterjacobian{counter}=funch{counter}{8}{7};
            OCMATINDIF.variationalobjectivefunctionderivativetime{counter}=funch{counter}{8}{8};
            OCMATINDIF.variationalsalvagevalue{counter}=funch{counter}{8}{9};
        end
        if ~isempty(OCMATINDIF.userbc) && OCMATINDIF.userbc{ii}
            OCMATINDIF.userfunctionbc{counter}=funch{counter}{5}{9};
        end

        if OCMATINDIF.exogenousfunction
            OCMATINDIF.exogenousdynamics{counter}=funch{counter}{1}{4};
            OCMATINDIF.exogenousjacobian{counter}=funch{counter}{2}{8};
            OCMATINDIF.exogenousparameterjacobian{counter}=funch{counter}{2}{9};
            OCMATINDIF.exogenousdct{counter}=funch{counter}{2}{10};
        end
        if OCMATFTE.variationalcalculation(counter)
            OCMATFTE.variationaldynamics{counter}=funch{counter}{4}{4};
            OCMATFTE.variationaljacobian{counter}=funch{counter}{4}{5};
            OCMATFTE.variationalparameterjacobian{counter}=funch{counter}{4}{6};
            OCMATFTE.variationalhamiltonian{counter}=funch{counter}{4}{10};
            OCMATFTE.variationalguard{counter}=funch{counter}{5}{9};
            OCMATFTE.variationalreset{counter}=funch{counter}{5}{10};
            OCMATFTE.variationalbcinitial{counter}=funch{counter}{5}{11};
            OCMATFTE.variationalbctransversality{counter}=funch{counter}{5}{12};
            if OCMATFTE.exogenousfunction
                OCMATFTE.exogenousvariationaldynamics{counter}=funch{counter}{4}{7};
                OCMATFTE.exogenousjacobian4variationalargument{counter}=funch{counter}{4}{8};
                OCMATFTE.exogenousvariationaldynamicsjacobian{counter}=funch{counter}{4}{9};
            end
        end
        %OCMATINDIF.canonicalsystemderivativetime{counter}=funch{counter}{2}{3};
        % general function
        OCMATINDIF.plotcontinuation{counter}=funch{counter}{11};
        OCMATINDIF.testadmissibility{counter}=funch{counter}{12};
        OCMATINDIF.datapath{counter}=funch{counter}{20};
        OCMATINDIF.saveintermediatefiles{counter}=funch{counter}{21};

        hybridinfo{counter}=OCMATINDIF.hybridinfo{counter}();

        for kk=1:numel(hybridinfo{counter}.arcarg)
            domaindata{counter}(kk)=OCMATINDIF.domain{counter}(hybridinfo{counter}.arcarg(kk));
        end
        for kk=1:numel(domaindata{counter})
            OCMATCONT.DOMAINDDATA{counter}(kk).numode=domaindata{counter}(kk).odedim;
            OCMATCONT.DOMAINDDATA{counter}(kk).numae=domaindata{counter}(kk).aedim;
            OCMATCONT.DOMAINDDATA{counter}(kk).daeorder=domaindata{counter}(kk).daeorder;
            OCMATCONT.DOMAINDDATA{counter}(kk).numeq=numel(domaindata{counter}(kk).daeorder);%number of equations
            if objectivevaluecalc
                OCMATCONT.DOMAINDDATA{counter}(kk).numode=domaindata{counter}(1).odedim+1;
                if jj==1
                    OCMATCONT.DOMAINDDATA{counter}(kk).numeq=OCMATCONT.DOMAINDDATA{counter}(1).numeq+1;%number of equations
                else
                    OCMATCONT.DOMAINDDATA{counter}(kk).numeq=OCMATCONT.DOMAINDDATA{counter}(1).numeq;%number of equations
                end
            end
            OCMATCONT.DOMAINDDATA{counter}(kk).eqcoord=1:OCMATCONT.DOMAINDDATA{counter}(kk).numeq;
            OCMATCONT.DOMAINDDATA{counter}(kk).odecoord=1:domaindata{counter}(kk).odedim;
            OCMATCONT.DOMAINDDATA{counter}(kk).aecoord=domaindata{counter}(kk).odedim+(1:domaindata{counter}(kk).aedim);
        end
    end
end

arccounter=0;
funccounter=0;
funcindex=[];
partindex=[];
relarcindex=[];
partcounter=0;
timeshift=[];
coordcounter=0;
for ii=1:OCMATINDIF.indifferenceorder
    arcn=[];
    arcarg=[];
    relpartcounter=0;
    relarccounter=0;
    if isoctrajectory(compTrj(ii))
        arcn{1}=arcnum(compTrj(ii));
        arcarg=arcargument(compTrj(ii));
    else
        arcn=arcnum(compTrj(ii));
        arcarg=arcargument(compTrj(ii));
        arcarg=[arcarg{:}];
    end
    arccounter0=arccounter+1;
    for jj=1:numberofparts(compTrj(ii))
        relpartcounter=relpartcounter+1;
        partcounter=partcounter+1;
        for kk=1:arcn{jj}
            relarccounter=relarccounter+1;
            arccounter=arccounter+1;
            coordcounter=coordcounter+1;
            if orderofmodel==1
                funcindex(arccounter)=jj;
            else
                funccounter=funccounter+1;
                funcindex(arccounter)=funccounter;
            end
            partindex(arccounter)=partcounter;
            relpartindex(arccounter)=relpartcounter;
            relarcindex(arccounter)=relarccounter;
            indifferenceindex(arccounter)=ii;

        end
        changepart(arccounter)=1;
    end
    changepart(arccounter)=-1;
    timeshift(ii+1)=arccounter;
    relarcargument{ii}=arcarg;
    relcoord{ii}=arccounter0:arccounter;
end
OCMATINDIF.indifferenceindex=indifferenceindex;
OCMATINDIF.funcindex=funcindex;
OCMATINDIF.partindex=partindex;
OCMATINDIF.relpartindex=relpartindex;
OCMATINDIF.relarcindex=relarcindex;
OCMATINDIF.totalnumberofparts=length(OCMATINDIF.partindex);
OCMATINDIF.totalnumberofarcs=arccounter;
OCMATINDIF.timeshift=timeshift(1:OCMATINDIF.indifferenceorder);
OCMATINDIF.changepart=changepart(1:end);
OCMATINDIF.relarcargument=relarcargument;
OCMATINDIF.relcoord=relcoord;

OCMATINDIF.freeparametervector=freeparametervector;
% counter=0;
% for ii=1:OCMATINDIF.indifferenceorder;
%     for jj=1:numberofparts(compTrj(ii))
%         counter=counter+1;
%         if orderofmodel==1
%             OCMATINDIF.parametervalue{ii}=parametervalue(ocObj);
%         else
%             if ismultistagemodel(ocObj(ii))
%                 OCMATINDIF.parametervalue{counter}=parametervalue(ocObj(ii));
%             else
%                 OCMATINDIF.parametervalue{counter}=parametervalue(ocObj(ii,jj));
%                 parametervaluetmp{ii}{jj}=parametervalue(ocObj(ii,jj));
%             end
%         end
%     end
% end
parametervaluetmp=parametervalue(ocObj);
counter=0;
for ii=1:OCMATINDIF.indifferenceorder
    for jj=1:nummodel(ii)
        counter=counter+1;
        OCMATINDIF.parametervalue{counter}=parametervaluetmp{ii}{jj};
    end
end
if ~isempty(freeparametervector)
    % generate freeparametervector as a cell array of size
    % {ii}{jj}, 1<=ii<=indifferenceorder, 1<=jj<=nummodel(ii)
    wrongcellsize=0;
    freeparametervectortmp=cell(1,OCMATINDIF.indifferenceorder);
    if iscell(freeparametervector) && numel(freeparametervector)==1
        freeparametervector=freeparametervector{1};
    end
    if iscell(freeparametervector)
        if numel(freeparametervector)==1
            if all(nummodel==numel(freeparametervector{1}))
                for ii=1:OCMATINDIF.indifferenceorder
                    freeparametervectortmp{ii}=freeparametervector;
                end
            else
                wrongcellsize=1;
            end
        elseif  numel(freeparametervector)==OCMATINDIF.indifferenceorder
            for ii=1:OCMATINDIF.indifferenceorder
                if nummodel(ii)==length(freeparametervector{ii})
                    freeparametervectortmp{ii}=freeparametervector{ii};
                else
                    wrongcellsize=1;
                end
            end
        else
            wrongcellsize=1;
        end
    elseif ischar(freeparametervector)
        % it is assumed that for each model the same parameters are free
        for ii=1:OCMATINDIF.indifferenceorder
            for jj=1:nummodel(ii)
                freeparametervectortmp{ii}{jj}=freeparametervector;
            end
        end
    else
        ocmatmsg('Free parameter vector has to be a cell-array or character.')
        return
    end
    if wrongcellsize
        ocmatmsg('Free parameter vector cell has wrong size.')
        return
    end
    freeparametervector=freeparametervectortmp;
    counter=0;
    for ii=1:OCMATINDIF.indifferenceorder
        ocObjIndiff=ocObj(ii);
        for jj=1:nummodel(ii)
            counter=counter+1;
            if ischar(freeparametervector{ii}{jj})
                OCMATINDIF.freeparameterindex{counter}=parameterindex(ocObjIndiff(jj),freeparametervector{ii}{jj});
                freeparametername{counter}=parametername(ocObjIndiff(jj),OCMATINDIF.freeparameterindex{counter});
                freeparametervalue{counter}=parametervalue(ocObjIndiff(jj),OCMATINDIF.freeparameterindex{counter});
            else
                OCMATINDIF.freeparameterindex{counter}=freeparametervector{ii}{jj};
                freeparametervalue{counter}=parametervalue(ocObjIndiff(jj),freeparametervector{ii}{jj});
            end
        end
    end
    numberoffreeparameter=-inf;
    freeparameternametotal='';
    freeparametervaluetotal=[];
    counter=0;
    for ii=1:OCMATINDIF.indifferenceorder
        for jj=1:nummodel(ii)
            counter=counter+1;
            numberoffreeparameter=max(numberoffreeparameter,length(OCMATINDIF.freeparameterindex{counter}));
            freeparameternametotal=[freeparameternametotal;cellstr(freeparametername{counter})];
            freeparametervaluetotal=[freeparametervaluetotal;freeparametervalue{counter}.'];
        end
    end
    %freeparametervalue=parametervaluetmp{1}{1}(OCMATINDIF.freeparameterindex{1});
    [freeparameternametotal,idx]=unique(freeparameternametotal);
    freeparametervalue=freeparametervaluetotal(idx);
    freeparametervalue=freeparametervalue(:).';
    counter=0;
    for ii=1:OCMATINDIF.indifferenceorder
        ocObjIndiff=ocObj(ii);
        for jj=1:nummodel(ii)
            counter=counter+1;
            actualnames=parametername(ocObjIndiff(jj),OCMATINDIF.freeparameterindex{counter});
            counter2=0;
            for kk=1:length(actualnames)
                idx=strmatch(actualnames{kk},freeparameternametotal);
                if ~isempty(idx)
                    counter2=counter2+1;
                    freeparametervectorcoordinate{counter}(counter2)=idx;
                end
            end
        end
    end
    %OCMATINDIF.parameterindex=totalparamterindex; %assumed that number and position of parameter values are the same for each model
    %     for ii=1:OCMATINDIF.indifferenceorder
    %         if orderofmodel==1
%             OCMATINDIF.parameterindex=parameterindex(ocObj,freeparametervector);
%             for jj=1:numberofparts(compTrj(ii))
%                 OCMATINDIF.freeinitialparametervalue{jj}=OCMATINDIF.parametervalue{jj}(OCMATINDIF.parameterindex(kk));
%             end
%         else
%             for jj=1:numberofparts(compTrj(ii))
%                 counter=counter+1;
%                 OCMATINDIF.freeinitialparametervalue{counter}=OCMATINDIF.parametervalue{counter}(OCMATINDIF.parameterindex);
%                 OCMATINDIF.freeinitialparametervalue{counter}=OCMATINDIF.freeinitialparametervalue{counter}(:);
%             end
%         end
%     end
end
% freestatevector is a matrix, where each column is a free vector and the
% number of rows is equal to the number of states
OCMATINDIF.freestatevector=freestatevector;



totcounter=0;
totarcint=[];
totarcarg=[];

for ii=1:OCMATINDIF.indifferenceorder
    arcint=arcinterval(compTrj(ii));
    arcarg=arcargument(compTrj(ii));
    arcn=[];
    if iscell(arcint)
        for jj=1:length(arcint)
            if jj==1
                arcinttmp=arcint{jj};
            else
                arcinttmp=[arcinttmp arcint{jj}(2:end)];
            end
        end
        arcint=arcinttmp;
        arcarg=[arcarg{:}];
        arcn=arcnum(compTrj(ii));

    else
        arcn{1}=arcnum(compTrj(ii));
    end
    totarcint=[totarcint arcint];
    totarcarg=[totarcarg arcarg];

    counterst=0;
    counter=1;
    counterct=0;
    counteret=0;
    counterstage=1;
    totcounter0=totcounter+1;
    switchtimeindex=[];
    conecttimeindex=[];
    timepointscharacterization=[];
    opttrans=optimaltransition(compTrj(ii));
    for jj=1:length(arcn)
        for kk=1:arcn{jj}
            counter=counter+1;
            if kk<length(arcn{jj})
                counterst=counterst+1;
                switchtimeindex(counterst)=counter;
                timepointscharacterization(counter)=1;
            else
                if jj<length(arcn)
                    counterct=counterct+1;
                    conecttimeindex(counterct)=counter;
                    if opttrans(counterstage)==1
                        timepointscharacterization(counter)=-1;
                    else
                        timepointscharacterization(counter)=0;
                    end
                else
                    counteret=counteret+1;
                    endtimeindex(counteret)=counter;
                    timepointscharacterization(counter)=0;
                end
                counterstage=counterstage+1;
            end
        end
        
    end
    totcounter=totcounter+counter;
    if isempty(connectiontimeindex)
        OCMATINDIF.connectiontimeindex{ii}=conecttimeindex;
    else
        OCMATINDIF.connectiontimeindex{ii}=connectiontimeindex{ii};
    end
    %OCMATINDIF.connectiontimeindex{ii}(opttrans(1:end-1)==1)=[];
    OCMATINDIF.opttimeindex{ii}=conecttimeindex(opttrans(1:end-1)==1);
    OCMATINDIF.switchtimeindex{ii}=switchtimeindex;
    OCMATINDIF.endtimeindex{ii}=endtimeindex;
    OCMATINDIF.arctimeindex{ii}=totcounter0:totcounter;
    OCMATINDIF.timepointscharacterization{ii}=timepointscharacterization; % 0 fixed time, 1 switching time, -1 optimized time, 2 continuation time
end
if isempty(freetimevector)
    OCMATINDIF.freetimevector=[];
end
OCMATINDIF.initialarcinterval=totarcint;
OCMATINDIF.totarcarg=totarcarg;
OCMATINDIF.initialstate=state(ocObj(1),compTrj(1),1);
OCMATINDIF.initialstate=OCMATINDIF.initialstate(:,1);
OCMATINDIF.fixinitialstate=fixinitialstate;
if ~isempty(fixinitialstate)
    OCMATINDIF.fixinitialstatecoordinate=fixinitialstate;
end
OCMATINDIF.fixendstate=fixendstate;
%OCMATINDIF.continuationindex=continuationindex;

OCMATINDIF.objectivevaluecalc=objectivevaluecalc;
OCMATINDIF.exogenousfunction=exogenousfunction;
if exogenousfunction
    OCMATINDIF.exogenousinitialstates=exogenousinitialstates;
    OCMATINDIF.exogenousnumberofstates=length(exogenousinitialstates);
end

switch continuationtype
    case 'initialstate'
        OCMATINDIF.continuationtype=0;
        OCMATINDIF.continuationvector=continuationtarget(:)-OCMATINDIF.initialstate(continuationindex);
        continuationparameter=0;
    case 'endstate'
        OCMATINDIF.continuationtype=0;
        OCMATINDIF.continuationvector=continuationtarget(:)-OCMATINDIF.endstate(continuationindex);
        continuationparameter=0;

    case 'time'
        OCMATINDIF.continuationtype=1;
        if length(continuationtarget)==1
            continuationtarget=repmat(continuationtarget,1,OCMATINDIF.indifferenceorder);
        end
        for ii=1:OCMATINDIF.indifferenceorder
            OCMATINDIF.continuationvector(ii)=continuationtarget(ii)-OCMATINDIF.initialarcinterval(OCMATINDIF.arctimeindex{ii}(continuationindex(ii)));
        end
        continuationparameter=0;
        OCMATINDIF.continuationindex=continuationindex;
    case 'parameter'
        % it is assumed that continuationtarget is of the same structure
        % as continuationindex
        OCMATINDIF.continuationtype=2;
        wrongcellsize=0;
        wrongcontinuationtargetstruct=0;
        continuationindextmp=cell(1,OCMATINDIF.indifferenceorder);
        continuationtargettmp=cell(1,OCMATINDIF.indifferenceorder);
        if iscell(continuationindex) && numel(continuationindex)==1
            if iscell(continuationtarget) && numel(continuationtarget)==1
            else
                wrongcontinuationtargetstruct=1;
            end
            continuationindex=continuationindex{1};
            continuationtarget=continuationtarget{1};
        end
        if iscell(continuationindex)
            if numel(continuationindex)==1
                if ~numel(continuationtarget)==1
                    wrongcontinuationtargetstruct=1;
                end
                if all(nummodel==numel(continuationindex{1}))
                    for ii=1:OCMATINDIF.indifferenceorder
                        continuationindextmp{ii}=continuationindex;
                        continuationtargettmp{ii}=continuationtarget;
                    end
                else
                    wrongcellsize=1;
                end
            elseif  numel(continuationindex)==OCMATINDIF.indifferenceorder
                for ii=1:OCMATINDIF.indifferenceorder
                    if nummodel(ii)==length(continuationindex{ii})
                        for jj=1:length(continuationindex{ii})
                            continuationindextmp{ii}{jj}=continuationindex{ii}{jj};
                            if isempty(continuationindex{ii}{jj})
                                continuationtargettmp{ii}{jj}=[];
                            else
                                if iscell(continuationtarget)
                                    continuationtargettmp{ii}{jj}=continuationtarget{ii}{jj};
                                else
                                    if length(continuationtarget)==length(continuationindex{ii})
                                        continuationtargettmp{ii}{jj}=continuationtarget(jj);
                                    else
                                        continuationtargettmp{ii}{jj}=continuationtarget;
                                    end
                                end
                            end
                        end
                    else
                        wrongcellsize=1;
                    end
                end
            else
                wrongcellsize=1;
            end
        elseif ischar(continuationindex)
            if ~length(continuationindex)==length(continuationtarget)
                    wrongcontinuationtargetstruct=1;
            end
            % it is assumed that for each model the same parameters are free
            for ii=1:OCMATINDIF.indifferenceorder
                for jj=1:nummodel(ii)
                    continuationindextmp{ii}{jj}=continuationindex;
                    continuationtargettmp{ii}{jj}=continuationtarget;
                end
            end
        else
            ocmatmsg('Parameter vector has to be a cell-array or character.')
            return
        end
        if wrongcellsize
            ocmatmsg('Parameter vector cell has wrong size.')
            return
        end
        if wrongcontinuationtargetstruct
            ocmatmsg('Target parameter vector cell has wrong structure.')
            return
        end
        continuationindex=continuationindextmp;
        continuationtarget=continuationtargettmp;
        counter=0;
        for ii=1:OCMATINDIF.indifferenceorder
            ocObjIndiff=ocObj(ii);
            for jj=1:nummodel(ii)
                counter=counter+1;
                if ischar(continuationindex{ii}{jj})
                    continuationindex{ii}{jj}=parameterindex(ocObjIndiff(jj),continuationindex{ii}{jj});
                    OCMATINDIF.continuationindex{counter}=continuationindex{ii}{jj};
                else
                    OCMATINDIF.continuationindex{counter}=continuationindex{ii}{jj};
                end
                if length(continuationtarget{ii}{jj})==length(continuationindex{ii}{jj})
                    OCMATINDIF.initialparameter{counter}=parametervaluetmp{ii}{jj}(continuationindex{ii}{jj});
                    if isempty(continuationtarget{ii}{jj})
                        OCMATINDIF.continuationvector{counter}=[];
                    else
                        OCMATINDIF.continuationvector{counter}=continuationtarget{ii}{jj}-OCMATINDIF.initialparameter{counter};
                    end
                else
                    ocmatmsg('Target parameter vector cell has wrong size.')
                    return
                end
            end
        end
        continuationparameter=0;
end
scoord=statecoord(ocObj(1));
scoord=scoord{1};
cscoord=costatecoord(ocObj(1));
cscoord=cscoord{1};

lastcoordinate=cscoord(end);
if OCMATINDIF.explicitconnectiontime
    OCMATINDIF.connectiontimenumber=connectiontimenumber(ocObj(1));
    OCMATINDIF.dhamiltoniandctcoord=lastcoordinate+(1:max(OCMATINDIF.connectiontimenumber)); %it is assumed that every part has the same number of states and costates
    OCMATINDIF.connectiontimenumber=max(OCMATINDIF.connectiontimenumber);

    for ii=1:OCMATINDIF.indifferenceorder
        ocTrjCell=[];
        dHdctval=zeros(OCMATINDIF.connectiontimenumber,1);
        ocTrjMP{ii}=compTrj(ii);
        for jj=1:numberofparts(compTrj(ii))
            ocTrjStruct=struct(ocTrjMP{ii}(jj));
            if length(ocTrjMP{ii}(jj).y(:,1))==lastcoordinate
                dSdctval=dsalvagedct(ocObj(jj),ocTrjMP);
                val=dhamiltoniandct(ocObj(jj),ocTrjMP,1);
                ocTrjStruct.y(OCMATINDIF.dhamiltoniandctcoord,:)=repmat(dSdctval,1,size(val,2))+repmat(dHdctval,1,size(val,2))+[zeros(OCMATINDIF.connectiontimenumber,1) cumsum(val(:,1:end-1)+val(:,2:end),1)/2.*repmat(diff(time(ocObj(jj),ocTrjMP,1)),OCMATINDIF.connectiontimenumber,1)];
                ocTrjCell{ii}{jj}=octrajectory(ocTrjStruct);
                dHdctval=ocTrjStruct.y(OCMATINDIF.dhamiltoniandctcoord,end);
            else
                ocTrjCell{jj}=ocTrjMP{ii}(jj);
                dHdctval=ocTrjMP{ii}(jj).y(OCMATINDIF.dhamiltoniandctcoord,end);
            end
        end
        ocTrjMP{ii}=mmultipath(ocTrjCell);
    end
    lastcoordinate=lastcoordinate+OCMATINDIF.connectiontimenumber;
end
compTrj=occomposite(ocTrjMP);
if objectivevaluecalc
    for ii=1:OCMATINDIF.indifferenceorder
        if ~isempty(objectivevaluecoordinate)
            actPath=compTrj(ii);
            for jj=1:numberofparts(compTrj(ii))
                objpath{jj}=actPath(jj).y(objectivevaluecoordinate,:);
            end
        else
            objpath=objectivevaluepath(compTrj(ii));
        end
        if isempty(objpath)
            % not correctly implemented now
            ctr=0;
            objpath=cell(1,numparts);
            for jj=1:OCMATINDIF.indifferenceorder
                Oval=0;
                for kk=1:nummod
                    ctr=ctr+1;
                    if numparts==nummod
                        actocObj=ocObj(ctr);
                    elseif numparts==nummod*OCMATINDIF.indifferenceorder
                        actocObj=ocObj(jj);
                    end
                    if length(compTrj(ctr).y(:,1))==2*statenum(actocObj)
                        o=objectivefunction(actocObj,compTrj(ctr),1);
                        objpath{ctr}=Oval+[0 cumsum((o(1:end-1)+o(2:end))/2.*diff(time(actocObj,compTrj(jj),1)))];
                        Oval=ocTrjStruct.y(end,end);
                    else
                        objpath{ctr}=compTrj(ctr).y(end,:);
                        Oval=objpath{ctr}(end);
                    end
                end
            end
            compTrj=addobjectivevaluepath(compTrj,objpath,lastcoordinate+1);
        end

    end
    OCMATINDIF.objectivevaluecoord=objectivevaluecoordinate;
end
if OCMATINDIF.exogenousfunction
    OCMATINDIF.exogenousdynamicscoordinate=OCMATINDIF.objectivevaluecoord+1:OCMATINDIF.objectivevaluecoord+OCMATINDIF.exogenousnumberofstates;
else
    OCMATINDIF.exogenousdynamicscoordinate=[];
end

OCMATINDIF.statecostatecoord=[scoord(:).' cscoord(:).'];
if objectivevaluecalc
    if OCMATINDIF.explicitconnectiontime
        OCMATINDIF.zeros2npctpobj=zeros(OCMATINDIF.statecostatecoord(end)+OCMATINDIF.connectiontimenumber,OCMATINDIF.connectiontimenumber);
    else
        OCMATINDIF.zeros2npctpobj=zeros(OCMATINDIF.statecostatecoord(end)+OCMATINDIF.connectiontimenumber,OCMATINDIF.connectiontimenumber);
    end
end
sol=generatesolstruct(compTrj);

if ~isempty(OCMATINDIF.freestatevector)
    OCMATINDIF.freestatevectorcoordinate=length(sol.parameters)+(1:size(OCMATINDIF.freestatevector,2));
    sol.parameters=[sol.parameters zeros(1,size(OCMATINDIF.freestatevector,2))];
end
if ~isempty(freeparametervector)
    OCMATINDIF.numberoffreeparameter=numberoffreeparameter;
    %OCMATINDIF.freeparametervectorcoordinate=length(sol.parameters)+(1:numberoffreeparameter);
    for ii=1:OCMATINDIF.totalnumberofparts
       OCMATINDIF.freeparametervectorcoordinate{ii}=length(sol.parameters)+freeparametervectorcoordinate{ii};
    end
    
    sol.parameters=[sol.parameters freeparametervalue];
end
for ii=1:OCMATINDIF.indifferenceorder
    OCMATINDIF.timepointsfreeindex{ii}=find(OCMATINDIF.timepointscharacterization{ii}==1);
    OCMATINDIF.timepointsfree4parcoordinate{ii}=length(sol.parameters)+(1:length(OCMATINDIF.timepointsfreeindex{ii}));
    sol.parameters=[sol.parameters sol.arcinterval(OCMATINDIF.timepointscharacterization{ii}==1)];

    OCMATINDIF.timepointsoptimizedindex{ii}=find(OCMATINDIF.timepointscharacterization{ii}==-1);
    OCMATINDIF.timepointsoptimize4parcoordinate{ii}=length(sol.parameters)+(1:length(OCMATINDIF.timepointsoptimizedindex{ii}));

    %OCMATINDIF.opttimecoordinate{ii}=length(sol.parameters)+(1:length(OCMATINDIF.opttimeindex{ii}));
    actarcinterval=sol.arcinterval(OCMATINDIF.arctimeindex{ii});
    sol.parameters=[sol.parameters actarcinterval(OCMATINDIF.timepointsoptimizedindex{ii})];
end
if ~isempty(freetimevector)
    if length(freetimevector)==1
        freetimevector=repmat(freetimevector,1,OCMATINDIF.indifferenceorder);
    end
    if  length(freetimevector)==OCMATINDIF.indifferenceorder
        OCMATINDIF.freetimevector=freetimevector;
    else
        return
    end
    OCMATINDIF.freetimevectorindex=length(sol.parameters)+1;
    freetime=zeros(1,OCMATINDIF.indifferenceorder);
    for ii=1:OCMATINDIF.indifferenceorder
        arcint=arcinterval(compTrj(ii));
        arcint=unique([arcint{:}]);
        freetime(ii)=arcint(freetimevector(ii));
    end
    freetime=unique(freetime);
    freetime=freetime(1);
    if length(freetime)==1
        sol.parameters=[sol.parameters freetime];
    else
        return
    end
end
OCMATINDIF.continuationcoordinate=length(sol.parameters)+1;
sol.parameters=[sol.parameters continuationparameter];

OCMATINDIF.statecoordinate=scoord;
OCMATINDIF.dFDO=[]; % derivative of the canonical system with respect to the objective variable
OCMATINDIF.dFDE=[]; % derivative of the canonical system with respect to the exogenous variable
OCMATINDIF.dFDCT=[]; % derivative of the canonical system with respect to the connecting times
OCMATINDIF.dFDPAR=[]; % derivative of the canonical system with respect to the connecting times
OCMATINDIF.dFODX=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODO=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODE=[]; % derivative of the objective dynamics with respect to the exogenous variable
OCMATINDIF.dFODCT=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFODPAR=[]; % derivative of the objective dynamics with respect to the objective variable
OCMATINDIF.dFEDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFEDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFEDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFEDCT=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFEDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFHCTDX=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFHCTDO=[]; % derivative of the exogenous dynamics with respect to the objective variable
OCMATINDIF.dFHCTDE=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFHCTDCT=[]; % derivative of the exogenous dynamics with respect to the exogenous variable
OCMATINDIF.dFHCTDPAR=[]; % derivative of the exogenous dynamics with respect to the exogenous variable

OCMATINDIF.dFDPAR=zeros(OCMATINDIF.statecostatecoord(end),length(sol.parameters));
if OCMATINDIF.explicitconnectiontime
    OCMATINDIF.dFDCT=zeros(OCMATINDIF.statecostatecoord(end),OCMATINDIF.connectiontimenumber);
    OCMATINDIF.dFHCTDX=zeros(OCMATINDIF.connectiontimenumber,OCMATINDIF.statecostatecoord(end));
    if objectivevaluecalc
        OCMATINDIF.dFODCT=zeros(1,OCMATINDIF.connectiontimenumber);
    end
    if exogenousfunction
        OCMATINDIF.dFEDCT=zeros(OCMATINDIF.exogenousnumberofstates,OCMATINDIF.connectiontimenumber);
    end
    OCMATINDIF.dFHCTDCT=zeros(OCMATINDIF.connectiontimenumber);
    OCMATINDIF.dFHCTDPAR=zeros(OCMATINDIF.connectiontimenumber,length(sol.parameters));
end

if objectivevaluecalc
    OCMATINDIF.dFDO=zeros(OCMATINDIF.statecostatecoord(end),1);
    OCMATINDIF.dFODO=0;
    OCMATINDIF.dFODX=zeros(1,OCMATINDIF.statecostatecoord(end));
    if OCMATINDIF.explicitconnectiontime
        OCMATINDIF.dFHCTDO=zeros(OCMATINDIF.connectiontimenumber,1);
    end
    if exogenousfunction
        OCMATINDIF.dFODE=zeros(1,OCMATINDIF.exogenousnumberofstates);
    end
    OCMATINDIF.dFODPAR=zeros(1,length(sol.parameters));
end
if OCMATINDIF.exogenousfunction
    OCMATINDIF.dFDE=zeros(OCMATINDIF.statecostatecoord(end),OCMATINDIF.exogenousnumberofstates);
    OCMATINDIF.dFEDE=zeros(OCMATINDIF.exogenousnumberofstates);
    OCMATINDIF.dFEDX=zeros(OCMATINDIF.exogenousnumberofstates,OCMATINDIF.statecostatecoord(end));
    if OCMATINDIF.objectivevaluecalc
        OCMATINDIF.dFEDO=zeros(OCMATINDIF.exogenousnumberofstates,1);
    end
    if OCMATINDIF.explicitconnectiontime
        OCMATINDIF.dFEDCT=zeros(OCMATINDIF.exogenousnumberofstates,OCMATINDIF.connectiontimenumber);
    end
    OCMATINDIF.dFEDPAR=zeros(OCMATINDIF.exogenousnumberofstates,length(sol.parameters));
end
OCMATINDIF.Jpar=[OCMATINDIF.dFDPAR;OCMATINDIF.dFHCTDPAR;OCMATINDIF.dFODPAR;OCMATINDIF.dFEDPAR];
 
OCMATINDIF.Jext=[[OCMATINDIF.dFDCT OCMATINDIF.dFDO OCMATINDIF.dFDE];[OCMATINDIF.dFHCTDCT OCMATINDIF.dFHCTDO OCMATINDIF.dFEDCT];[OCMATINDIF.dFODCT OCMATINDIF.dFODO OCMATINDIF.dFODE];[OCMATINDIF.dFEDCT OCMATINDIF.dFEDO OCMATINDIF.dFEDE]];
pathname=OCMATINDIF.datapath{1}();
[resultfile,globalvarfile]=OCMATINDIF.saveintermediatefiles{1}();
OCMATINDIF.basicresultfilename=fullocmatfile(pathname,resultfile);
OCMATINDIF.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCMATCONT.HE.numinitialcondition=[];

OCMATCONT.codimension=1;



function sol=generatesolstruct(ocMultiPath)
global OCMATINDIF
sol.x=[];
sol.y=[];
sol.arcarg=[];
sol.arcinterval=[];
ctr=0;
for ii=1:OCMATINDIF.indifferenceorder
    ocTrj=ocMultiPath(ii);
    for jj=1:numberofparts(ocTrj)
        sol.x=[sol.x independentvar(ocTrj(jj))+ctr];
        sol.y=[sol.y dependentvar(ocTrj(jj))];
        sol.arcarg=[sol.arcarg arcargument(ocTrj(jj))];
        arcint=arcinterval(ocTrj(jj));
        if jj==1
            sol.arcinterval=[sol.arcinterval arcint];
        else
            sol.arcinterval=[sol.arcinterval arcint(2:end)];
        end
        ctr=ctr+1;
    end
end
sol.parameters=[];
sol.x0=sol.arcinterval(1);
sol.solver='';
sol.parameters=[];
sol.solverinfo.tangent=[];
sol.solverinfo.coeff=[];
sol.solverinfo.tmesh=[];
