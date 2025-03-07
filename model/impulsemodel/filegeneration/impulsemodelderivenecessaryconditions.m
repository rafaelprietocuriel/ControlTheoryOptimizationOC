function ocStruct=impulsemodelderivenecessaryconditions(ocStruct,symkernel,opt)
%
% derive the necessary optimality conditions from the initialization file
if nargin==2
    opt=defaultocoptions;
end
JacobianCalc=getocoptions(opt,'INIT','Jacobian');
HessianCalc=getocoptions(opt,'INIT','Hessian');
if ~isstruct(ocStruct)
    ocStruct=loadmodeldata(modelname(ocStruct));
end
% general properties
arcidentifier=retrieveimpulsemodelinformation(ocStruct,'arcidentifier');
inequalitycontrolconstraintnum=retrieveimpulsemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
arcnum=retrieveimpulsemodelinformation(ocStruct,'arcnum');

%
% The Pontryagin function and its derivatives
%
ocStruct.pontryaginfunction.identifier='PF';
P=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunction');
Px=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDx',[],symkernel);
Pu=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDu',[],symkernel);
Px2=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDx2',[],symkernel);
Pu2=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDu2',[],symkernel);
Pux=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDuDx',[],symkernel);
ocStruct.pontryaginfunction.term=P.value;
ocStruct.pontryaginfunction.arcidentifier=arcidentifier.value;
ocStruct.pontryaginfunction.derivative.Dx.term=Px.value;
ocStruct.pontryaginfunction.derivative.Dx.identifier='DPFDx';
ocStruct.pontryaginfunction.derivative.Du.term=Pu.value;
ocStruct.pontryaginfunction.derivative.Du.identifier='DPFDu';
ocStruct.pontryaginfunction.derivative.Dx2.term=Px2.value;
ocStruct.pontryaginfunction.derivative.Dx2.identifier='DPFDx2';
ocStruct.pontryaginfunction.derivative.Du2.term=Pu2.value;
ocStruct.pontryaginfunction.derivative.Du2.identifier='DPFDu2';
ocStruct.pontryaginfunction.derivative.DuDx.term=Pux.value;
ocStruct.pontryaginfunction.derivative.DuDx.identifier='DPFDuDx';

%
% The Impulse Pontryagin function and its derivatives
%
ocStruct.impulsepontryaginfunction.identifier='IPF';
IP=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunction');
IPx=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunctionDx',[],symkernel);
IPu=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunctionDu',[],symkernel);
IPtau=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunctionDtau',[],symkernel);
% IPx2=retrieveimpulsemodelinformation(ocStruct,'impulsepontryaginfunctionDx2',[],symkernel);
% IPu2=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDu2',[],symkernel);
% IPux=retrieveimpulsemodelinformation(ocStruct,'pontryaginfunctionDuDx',[],symkernel);
ocStruct.impulsepontryaginfunction.term=IP.value;
ocStruct.impulsepontryaginfunction.arcidentifier=arcidentifier.value;
ocStruct.impulsepontryaginfunction.derivative.Dx.term=IPx.value;
ocStruct.impulsepontryaginfunction.derivative.Dx.identifier='DIPFDx';
ocStruct.impulsepontryaginfunction.derivative.Du.term=IPu.value;
ocStruct.impulsepontryaginfunction.derivative.Du.identifier='DIPFDu';
ocStruct.impulsepontryaginfunction.derivative.Dtau.term=IPtau.value;
ocStruct.impulsepontryaginfunction.derivative.Dtau.identifier='DIPFDtau';
% ocStruct.impulsepontryaginfunction.derivative.Dx2.term=IPx2.value;
% ocStruct.impulsepontryaginfunction.derivative.Dx2.identifier='DIPFDx2';
% ocStruct.impulsepontryaginfunction.derivative.Du2.term=IPu2.value;
% ocStruct.impulsepontryaginfunction.derivative.Du2.identifier='DIPFDu2';
% ocStruct.impulsepontryaginfunction.derivative.DuDx.term=IPux.value;
% ocStruct.impulsepontryaginfunction.derivative.DuDx.identifier='DIPFDuDx';

%
% The adjoint system
%
% determine possible algebraic equations
dldt=retrieveimpulsemodelinformation(ocStruct,'adjointsystem',arcidentifier.value,symkernel);
arcid=arcidentifier.value;
for ii=arcnum.value:-1:1
    ocStruct.foc.adjointsystem.algebraicequation.(arcidentifier2field(arcid{ii})).term=[];
end
ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).ode.term=dldt.value;
ocStruct.foc.adjointsystem.algebraicequation=orderfields(ocStruct.foc.adjointsystem.algebraicequation);

adjointevent=retrieveimpulsemodelinformation(ocStruct,'adjointevent',arcidentifier.value,symkernel);
ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).evt.term=adjointevent.value;
%
%
% Hamiltonian maximizing condition
%
% solve the corresponding first order necessary conditions and let user
% decide which solutions should be used if multiple solutions occur
%
arcwitheequalpontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'arcwitheequalpontryaginfunction');
maximizingvariable=retrieveimpulsemodelinformation(ocStruct,'maximizingvariable');
solutionindex=retrieveimpulsemodelinformation(ocStruct,'solutionindex');
for ii=1:numel(arcwitheequalpontryaginfunction.value)
    arcid=arcwitheequalpontryaginfunction.value{ii}{1};
    sol=retrieveimpulsemodelinformation(ocStruct,'maximizingsolution',arcid,symkernel);
    if numel(sol.value)>1
        solutionindices=cell2vectorstring(solutionindex.value{ii});
        % multiple solutions exist
        fprintf('Arc identifiers : %s\n\n',arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}))
        for jj=1:numel(sol.value)
            fprintf('Solution Number %d\n',jj)
            for kk=1:numel(maximizingvariable.value)
                fprintf('\t%s = %s\n',maximizingvariable.value{kk},sol.value(jj).(maximizingvariable.value{kk}))
            end
            fprintf('\n\n')
        end
        fprintf('%d solutions can be kept?\nSolutions chosen by default are %s\n ',numel(solutionindex.value{ii}),solutionindices)
        while 1
            accept=lower(input('Accept: (y)/n : ','s'));
            if isempty(accept)
                accept='y';
            end
            if strcmp(accept,'y')
                sol.value=eval(['sol.value(' solutionindices ')']);
                break
            elseif strcmp(accept,'n')
                break
            elseif regexp(accept,'q')
                break
            end
        end
        counter=0;
        while 1 && ~strcmp(accept,'y')
            counter=counter+1;
            fprintf('\n\n')
            answer=input('Input solution indices in Matlab notation : ','s');
            solutionindices=eval(['[' answer ']']);
            try
                if numel(solutionindices)~=numel(solutionindex.value{ii})
                    fprintf('Wrong number of solutions. Try again.\n')
                elseif ~all(ismember(solutionindices,1:numel(sol.value)))
                    fprintf('Solution indices are out of solution range.\n')
                else
                    sol.value=sol.value(solutionindices);
                    break
                end
            catch
                fprintf('\nWrong format %s. Try again.\n',answer)
            end
            if counter>2
                fprintf('\nGive up after %d try.\n',counter)
            end
        end
    end
    if ~isempty(sol.value)
        for jj=1:numel(arcwitheequalpontryaginfunction.value{ii})
            ocStruct.foc.value.control.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).control;
            if inequalitycontrolconstraintnum.value
                ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).lagrangemultcc;
            end
        end
    else
        for jj=1:numel(arcwitheequalpontryaginfunction.value{ii})
            ocStruct.foc.value.control.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term='';
        end
    end
end

arcwitheequalimpulsepontryaginfunction=retrieveimpulsemodelinformation(ocStruct,'arcwitheequalimpulsepontryaginfunction');
for ii=1:numel(arcwitheequalimpulsepontryaginfunction.value)
    arcid=arcwitheequalimpulsepontryaginfunction.value{ii}{1};
    sol=retrieveimpulsemodelinformation(ocStruct,'impulsemaximizingsolution','0',symkernel);
    if ~isempty(sol.value)
        for jj=1:numel(arcwitheequalimpulsepontryaginfunction.value{ii})
            ocStruct.foc.value.icontrol.(arcidentifier2field(arcwitheequalimpulsepontryaginfunction.value{ii}{jj})).term=sol.value(jj).icontrol;
        end
    end
end
% sort arc identifiers
arcargstring=char(fieldnames(ocStruct.foc.value.control));
[dum idx]=sort(str2num(arcargstring(:,4:end)));
ocStruct.foc.value.control=orderfields(ocStruct.foc.value.control,idx);
ocStruct.foc.value.icontrol=orderfields(ocStruct.foc.value.icontrol,idx);
if inequalitycontrolconstraintnum.value
    ocStruct.foc.value.lagrangemultcc=orderfields(ocStruct.foc.value.lagrangemultcc,idx);
end

% Transversality condition
DsalvagevalueDx=retrieveimpulsemodelinformation(ocStruct,'DsalvagevalueDx',[],symkernel);
DsalvagevalueDT=retrieveimpulsemodelinformation(ocStruct,'DsalvagevalueDT',[],symkernel);
ocStruct.foc.adjointsystem.DSDX=DsalvagevalueDx.value;
ocStruct.foc.adjointsystem.DSDT=DsalvagevalueDT.value;

%
% Jacobian of the canonical system
%
for ii=1:arcnum.value
    arc=arcidentifier2field(arcidentifier.value{ii});
    if strcmpi(JacobianCalc,'explicit')
        statename=retrieveimpulsemodelinformation(ocStruct,'statename',arcidentifier.value{ii});
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename',arcidentifier.value{ii});
        dxdt=retrieveimpulsemodelinformation(ocStruct,'specificstatedynamics',arcidentifier.value{ii},symkernel);
        dldt=retrieveimpulsemodelinformation(ocStruct,'specificadjointsystem',arcidentifier.value{ii},symkernel);
        dXdt=cell2vectorstring(term2optterm_impulsemodel(ocStruct,[dxdt.value dldt.value],arcidentifier.value{ii},opt));
        variable=cell2vectorstring([statename.value costatename.value]);
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');

        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername',arcidentifier.value{ii});
        variable=cell2vectorstring(parametername.value);
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');
    elseif strcmpi(JacobianCalc,'numerical')
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term='';
    else
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term='';
    end
end

%
% Jacobian of the events
%
for ii=1:arcnum.value
    arc=arcidentifier2field(arcidentifier.value{ii});
    if strcmpi(JacobianCalc,'explicit')
        statename=retrieveimpulsemodelinformation(ocStruct,'statename',arcidentifier.value{ii});
        costatename=retrieveimpulsemodelinformation(ocStruct,'costatename',arcidentifier.value{ii});
        stateevent=retrieveimpulsemodelinformation(ocStruct,'stateevent',arcidentifier.value{ii},symkernel);
        costateevent=retrieveimpulsemodelinformation(ocStruct,'adjointevent',arcidentifier.value{ii},symkernel);
        Xevent=cell2vectorstring(term2optterm_impulsemodel(ocStruct,[stateevent.value costateevent.value],arcidentifier.value{ii},opt));
        variable=cell2vectorstring([statename.value costatename.value]);
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDX.term=string2cell(removematrixstring(ocmatjacobian(Xevent,variable,symkernel,opt)),'matrix');

        parametername=retrieveimpulsemodelinformation(ocStruct,'parametername',arcidentifier.value{ii});
        variable=cell2vectorstring(parametername.value);
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDP.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');
        impulsetime=retrieveimpulsemodelinformation(ocStruct,'impulsetime');
        variable=cell2vectorstring(impulsetime.value);
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDtau.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');
    elseif strcmpi(JacobianCalc,'numerical')
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDP.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDtau.term='';
    else
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDP.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DEDtau.term='';
    end
end
%
% Hessian of the canonical system
%
% for ii=1:arcnum.value
%     if strcmpi(HessianCalc,'explicit')
%         arc=arcidentifier2field(arcidentifier.value{ii});
%         H=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier.value{ii},symkernel);
%         ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term=H.value;
%         Jxp=retrieveimpulsemodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier.value{ii},symkernel);
%         ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term=Jxp.value;
%     elseif strcmpi(HessianCalc,'numerical')
%         ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term='';
%         ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term='';
%     else
%         ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term='';
%         ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term='';
%     end
% end
