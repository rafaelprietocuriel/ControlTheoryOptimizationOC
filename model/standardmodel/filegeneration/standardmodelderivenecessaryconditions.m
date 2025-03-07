function ocStruct=standardmodelderivenecessaryconditions(ocStruct,symkernel,opt)
%
% derive the necessary optimality conditions from the initialization file
if nargin==2
    opt=defaultocoptions;
end
JacobianCalc=getocoptions(opt,'INIT','Jacobian');
HessianCalc=getocoptions(opt,'INIT','Hessian');
if implicitcontrols(ocStruct) && (strcmpi(JacobianCalc,'explicit') || strcmpi(HessianCalc,'explicit'))
    fprintf('\nThe symbolic derivation of the Jacobian/Hessian in case of implicit controls can take a considerably long time.\n')
    fprintf('To calculate the Jacobian/Hessian numerically set\n\n')
    fprintf('opt=setocoptions(''INIT'',''Jacobian'',''numerical'',''Hessian'',''numerical'');\n\n')
    fprintf('and restart the initial processing:\n\n')
    fprintf('ocStruct=processinitfile(''%s'',opt);\n\n',ocStruct.modelname)
    answer=input('Interrupt the processing (y)/n: ','s');
    while 1
        if isempty(answer)
            % default value 'y'
            answer='y';
        end
        if strcmpi(answer,'n')
            break
        elseif strcmpi(answer,'y')
            return
        end
    end
end
if ~isstruct(ocStruct)
    ocStruct=loadmodeldata(modelname(ocStruct));
end
% general properties
arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
%switchingnum=retrievemodelinformation(ocStruct,'switchingnum');
if inequalitystateconstraintnum.value
    statename=retrievemodelinformation(ocStruct,'statename');
    statenum=retrievemodelinformation(ocStruct,'statenum');
    statedynamics=retrievemodelinformation(ocStruct,'statedynamics');
    inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
    inequalitystateconstraintorder=retrievemodelinformation(ocStruct,'inequalitystateconstraintorder');
    if length(inequalitystateconstraintorder.value)~=inequalitystateconstraintnum.value
        ocmaterror('Entry for order of state constraint is missing.')
    end
    statevariable=statename.value;
    for ii=1:inequalitystateconstraintnum.value
        dhdt{ii}{1}=inequalitystateconstraint.value{ii};
        for jj=1:inequalitystateconstraintorder.value(ii)
            dhdt{ii}{jj+1}=0;
            for kk=1:statenum.value
                dhdt{ii}{jj+1}=dhdt{ii}{jj+1}+mystr2sym(ocmatdiff(dhdt{ii}{jj},statevariable{kk},symkernel))*mystr2sym(statedynamics.value{kk});
            end
            dhdt{ii}{jj+1}=char(dhdt{ii}{jj+1});
        end
        dhdt{ii}(1)=[];
    end
    ocStruct.constraint.function.state.timederivative=dhdt;
end
arcnum=retrievemodelinformation(ocStruct,'arcnum');

%
% The Pontryagin function and its derivatives
%
ocStruct.pontryaginfunction.identifier='PF';
P=retrievemodelinformation(ocStruct,'pontryaginfunction');
Px=retrievemodelinformation(ocStruct,'pontryaginfunctionDx',[],symkernel);
Pu=retrievemodelinformation(ocStruct,'pontryaginfunctionDu',[],symkernel);
Px2=retrievemodelinformation(ocStruct,'pontryaginfunctionDx2',[],symkernel);
Pu2=retrievemodelinformation(ocStruct,'pontryaginfunctionDu2',[],symkernel);
Pux=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDx',[],symkernel);
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
if implicitcontrols(ocStruct)
    Pulmmc=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDlmmc',[],symkernel);
    ocStruct.pontryaginfunction.derivative.DuDlmmc.term=Pulmmc.value;
    ocStruct.pontryaginfunction.derivative.DuDlmmc.identifier='DPFDuDlmmc';
    PuX=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDX',[],symkernel);
    ocStruct.pontryaginfunction.derivative.DuDX.term=PuX.value;
    ocStruct.pontryaginfunction.derivative.DuDX.identifier='DPFDuDX';
end

%
% The adjoint system
%
% determine possible algebraic equations
dldt=retrievemodelinformation(ocStruct,'adjointsystem',arcidentifier.value,symkernel);
arcid=arcidentifier.value;
for ii=arcnum.value:-1:1
    ae=retrievemodelinformation(ocStruct,'algebraicequation',arcid{ii},symkernel);
    ocStruct.foc.adjointsystem.algebraicequation.(arcidentifier2field(arcid{ii})).term=ae.value;
end
ocStruct.foc.adjointsystem.dynamics.(arcidentifier2field(arcidentifier.value)).ode.term=dldt.value;
%ocStruct.foc.adjointsystem.algebraicequation.(arcidentifier2field(arcid)).term=[];
ocStruct.foc.adjointsystem.algebraicequation=orderfields(ocStruct.foc.adjointsystem.algebraicequation);

%
%
% Hamiltonian maximizing condition
%
% solve the corresponding first order necessary conditions and let user
% decide which solutions should be used if multiple solutions occur
%
arcwitheequalpontryaginfunction=retrievemodelinformation(ocStruct,'arcwitheequalpontryaginfunction');
maximizingvariable=retrievemodelinformation(ocStruct,'maximizingvariable');
solutionindex=retrievemodelinformation(ocStruct,'solutionindex');
if ~isempty(maximizingvariable.value)
    for ii=1:numel(arcwitheequalpontryaginfunction.value)
        arcid=arcwitheequalpontryaginfunction.value{ii}{1};
        sol=retrievemodelinformation(ocStruct,'maximizingsolution',arcid,symkernel);
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
        for jj=1:numel(arcwitheequalpontryaginfunction.value{ii})
            ocStruct.foc.value.control.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).control;
            if inequalitycontrolconstraintnum.value
                ocStruct.foc.value.lagrangemultcc.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).lagrangemultcc;
            end
            if inequalitystateconstraintnum.value
                explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',arcwitheequalpontryaginfunction.value{ii}{jj},symkernel);
                ocStruct.foc.value.lagrangemultsc.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=sol.value(jj).lagrangemultsc;
                ocStruct.foc.value.state.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=explicitstatevalue.value;
            end
        end
    end
else
    for ii=1:numel(arcwitheequalpontryaginfunction.value)
        for jj=1:numel(arcwitheequalpontryaginfunction.value{ii})
            ocStruct.foc.value.control.(arcidentifier2field(arcwitheequalpontryaginfunction.value{ii}{jj})).term=[];
        end
    end
end
% sort arc identifiers
arcargstring=char(fieldnames(ocStruct.foc.value.control));
[dum,idx]=sort(str2num(arcargstring(:,4:end)));
ocStruct.foc.value.control=orderfields(ocStruct.foc.value.control,idx);
if inequalitycontrolconstraintnum.value
    ocStruct.foc.value.lagrangemultcc=orderfields(ocStruct.foc.value.lagrangemultcc,idx);
end
if inequalitystateconstraintnum.value
    ocStruct.foc.value.lagrangemultsc=orderfields(ocStruct.foc.value.lagrangemultsc,idx);
end

%
% Jacobian of the canonical system
%
for ii=1:arcnum.value
    arc=arcidentifier2field(arcidentifier.value{ii});
    if strcmpi(JacobianCalc,'explicit')
        if ~implicitcontrols(ocStruct)
            statename=retrievemodelinformation(ocStruct,'statename',arcidentifier.value{ii});
            costatename=retrievemodelinformation(ocStruct,'costatename',arcidentifier.value{ii});
            dxdt=retrievemodelinformation(ocStruct,'specificstatedynamics',arcidentifier.value{ii},symkernel);
            dldt=retrievemodelinformation(ocStruct,'specificadjointsystem',arcidentifier.value{ii},symkernel);
            dXdt=cell2vectorstring(term2optterm_standardmodel(ocStruct,[dxdt.value dldt.value],arcidentifier.value{ii},opt));
            variable=cell2vectorstring([statename.value costatename.value]);
            ocStruct.foc.canonicalsystem.derivative.(arc).DX.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');

            parametername=retrievemodelinformation(ocStruct,'parametername',arcidentifier.value{ii});
            variable=cell2vectorstring(parametername.value);
            ocStruct.foc.canonicalsystem.derivative.(arc).DP.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');
        else
            implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            implicitnonlinearcontrol=retrievemodelinformation(ocStruct,'implicitnonlinearcontrol',arcidentifier.value{ii});
            statename=retrievemodelinformation(ocStruct,'statename',arcidentifier.value{ii});
            costatename=retrievemodelinformation(ocStruct,'costatename',arcidentifier.value{ii});
            dxdt=retrievemodelinformation(ocStruct,'specificstatedynamics',arcidentifier.value{ii},symkernel);
            dldt=retrievemodelinformation(ocStruct,'specificadjointsystem',arcidentifier.value{ii},symkernel);
            dudt=retrievemodelinformation(ocStruct,'fulloptimalcontroldynamics',arcidentifier.value{ii},symkernel);
            dudt.value=dudt.value(implicitnonlinearcontrolindex.value);
            dXdt=cell2vectorstring(term2optterm_standardmodel(ocStruct,[dxdt.value dldt.value dudt.value],arcidentifier.value{ii},opt));
            variable=cell2vectorstring([statename.value costatename.value implicitnonlinearcontrol.value]);
            ocStruct.foc.canonicalsystem.derivative.(arc).DX.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');

            parametername=retrievemodelinformation(ocStruct,'parametername',arcidentifier.value{ii});
            variable=cell2vectorstring(parametername.value);
            ocStruct.foc.canonicalsystem.derivative.(arc).DP.term=string2cell(removematrixstring(ocmatjacobian(dXdt,variable,symkernel,opt)),'matrix');

        end
    elseif strcmpi(JacobianCalc,'numerical')
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term='';
    else
        ocStruct.foc.canonicalsystem.derivative.(arc).DX.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DP.term='';
    end
end

%
% Hessian of the canonical system
%
for ii=1:arcnum.value
    if strcmpi(HessianCalc,'explicit')
        arc=arcidentifier2field(arcidentifier.value{ii});
        H=retrievemodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier.value{ii},symkernel);
        ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term=H.value;
        Jxp=retrievemodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier.value{ii},symkernel);
        ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term=Jxp.value;
    elseif strcmpi(HessianCalc,'numerical')
        ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term='';
    else
        ocStruct.foc.canonicalsystem.derivative.(arc).DX2.term='';
        ocStruct.foc.canonicalsystem.derivative.(arc).DXP.term='';
    end
end
