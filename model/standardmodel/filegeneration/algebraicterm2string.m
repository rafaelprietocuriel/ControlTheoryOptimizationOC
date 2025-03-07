function algebraic=algebraicterm2string(ocStruct,classflag,vectflag,varargin)
%
% ALGEBRAICTERM2STRING process the algebraic terms of the structure
% ocStruct describing the model into strings which are printed into the
% model files


algebraic.term='';
algebraic.type='';
algebraic.info=[];
arcfield='';
classtype='';

if isempty(ocStruct)
    return
end

if nargin==2
    vectflag=0;
end

if nargin>=4
    arcfield=varargin{1};
end

if nargin>=5
    classtype=varargin{2};
end

switch classflag
    case 'controlvalue'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        if ~controlnum.value
            algebraic.term{1}='[]';
        else
            for ii=1:controlnum.value
                if ii==1
                    if controlnum.value>1
                        algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii}) '; ...'];
                    else
                        algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii})];
                    end
                elseif ii<controlnum.value
                    algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii})];
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'optimalvalue4statecontrol'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue4statecontrol',field2arcidentifier(arcfield));
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=['[' optimalvalue.value.(controlname.value{ii})];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii}) '; ...'];
            else
                algebraic.term{ii}=[optimalvalue.value.(controlname.value{ii})];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        
    case 'optimalcostatevalue'
        optimalcostatevalue=retrievemodelinformation(ocStruct,'optimalcostatevalue',field2arcidentifier(arcfield),getsymkernel);
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        for ii=1:costatenum.value
            if ii==1
                if costatenum.value>1
                    algebraic.term{ii}=['[' optimalcostatevalue.value.(costatename.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=optimalcostatevalue.value.(costatename.value{ii});
                end
            elseif ii<costatenum.value
                algebraic.term{ii}=[optimalcostatevalue.value.(costatename.value{ii}) '; ...'];
            else
                algebraic.term{ii}=[optimalcostatevalue.value.(costatename.value{ii})];
            end
        end
        if costatenum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
        
    case 'explicitstatevalue'
        explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',field2arcidentifier(arcfield),getsymkernel);
        statenum=retrievemodelinformation(ocStruct,'statenum');
        for ii=1:statenum.value
            if ii==1
                if statenum.value>1
                    algebraic.term{ii}=['[' explicitstatevalue.value{ii} '; ...'];
                else
                    algebraic.term{ii}=explicitstatevalue.value{ii};
                end
            elseif ii<statenum.value
                algebraic.term{ii}=[explicitstatevalue.value{ii} '; ...'];
            else
                algebraic.term{ii}=[explicitstatevalue.value{ii}];
            end
        end
        if statenum.value>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
        
    case 'symbolicstatevalue'
        explicitstatevalue=retrievemodelinformation(ocStruct,'explicitstatevalue',field2arcidentifier(arcfield),getsymkernel);
        statenum=retrievemodelinformation(ocStruct,'statenum');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:statenum.value
            if ii==1
                if statenum.value>1
                    algebraic.term{ii}=[symstr '([''[' explicitstatevalue.value{ii} ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' explicitstatevalue.value{ii} ']''])'];
                end
            elseif ii<statenum.value
                algebraic.term{ii}=['''' explicitstatevalue.value{ii} ','' ...'];
            else
                algebraic.term{ii}=['''' explicitstatevalue.value{ii} ']''])'];
            end
        end
        
    case 'equilibriumequation'
        equilibriumequation=retrievemodelinformation(ocStruct,'equilibriumequation',field2arcidentifier(arcfield),getsymkernel);
        equationnum=length(equilibriumequation.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:equationnum
            if ii==1
                if equationnum>1
                    algebraic.term{ii}=[symstr '([''[' equilibriumequation.value{ii} ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' equilibriumequation.value{ii} ']''])'];
                end
            elseif ii<equationnum
                algebraic.term{ii}=['''' equilibriumequation.value{ii} ','' ...'];
            else
                algebraic.term{ii}=['''' equilibriumequation.value{ii} ']''])'];
            end
        end
        
    case 'symboliccontrolvalue'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(controlname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(controlname.value{ii}) ']''])'];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ']''])'];
            end
        end
        
    case 'symboliccontrolvalue4statecontrol'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue4statecontrol',field2arcidentifier(arcfield));
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
         if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
       for ii=1:controlnum.value
            if ii==1
                if controlnum.value>1
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(controlname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(controlname.value{ii}) ']''])'];
                end
            elseif ii<controlnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(controlname.value{ii}) ']''])'];
            end
        end
    case 'symboliccostatevalue'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalcostatevalue',field2arcidentifier(arcfield),getsymkernel);
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:costatenum.value
            if ii==1
                if costatenum.value>1
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(costatename.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(costatename.value{ii}) ']''])'];
                end
            elseif ii<costatenum.value
                algebraic.term{ii}=['''' optimalvalue.value.(costatename.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(costatename.value{ii}) ']''])'];
            end
        end

    case 'lagrangemultcc'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        for ii=1:inequalitycontrolconstraintnum.value
            if ii==1
                if inequalitycontrolconstraintnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii})];
                end
            elseif ii<inequalitycontrolconstraintnum.value
                algebraic.term{ii}=[optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) '; ...'];
            else
                algebraic.term{ii}=optimalvalue.value.(lagrangemultipliercontrolname.value{ii});
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'symboliclagrangemultcc'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:inequalitycontrolconstraintnum.value
            if ii==1
                if inequalitycontrolconstraintnum.value>1
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ']''])'];
                end

            elseif ii<inequalitycontrolconstraintnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultipliercontrolname.value{ii}) ']''])'];
            end
        end

    case 'lagrangemultsc'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        lagrangemultscdt=retrievemodelinformation(ocStruct,'lagrangemultscdt',field2arcidentifier(arcfield),getsymkernel);
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        ctr=0;
        for ii=1:inequalitystateconstraintnum.value
            ctr=ctr+1;
            if ii==1
                algebraic.term{ctr}=['[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) '; ...'];
            else
                algebraic.term{ctr}=[optimalvalue.value.(lagrangemultiplierstatename.value{ii}) '; ...'];
            end
        end
        for ii=1:inequalitystateconstraintnum.value
            ctr=ctr+1;
            if ii<inequalitystateconstraintnum.value
                algebraic.term{ctr}=[lagrangemultscdt.value{ii} '; ...'];
            else
                algebraic.term{ctr}=lagrangemultscdt.value{ii};
            end
        end
        algebraic.term{ctr}=[algebraic.term{ctr} ']'];

    case 'symboliclagrangemultsc'
        optimalvalue=retrievemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrievemodelinformation(ocStruct,'lagrangemultiplierstatename');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:inequalitystateconstraintnum.value
            if ii==1
                if inequalitystateconstraintnum.value>1
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ','' ...'];
                else
                    algebraic.term{ii}=[symstr '([''[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ']''])'];
                end

            elseif ii<inequalitystateconstraintnum.value
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ','' ...'];
            else
                algebraic.term{ii}=['''' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) ']''])'];
            end
        end

    case 'objective'
        objectiveintegrand=retrievemodelinformation(ocStruct,'objectiveintegrand');
        switch classtype
            case 'current'
                discountfactor=retrievemodelinformation(ocStruct,'discountfactor');
                algebraic.term=[discountfactor.value '*(' objectiveintegrand.value ')'];
            case 'present'
                algebraic.term=objectiveintegrand.value;
        end

    case 'salvagevalue'
        salvagevalue=retrievemodelinformation(ocStruct,'salvagevalue');
        algebraic.term=salvagevalue.value;
        %endtimediscountfactor=retrievemodelinformation(ocStruct,'endtimediscountfactor');
        %algebraic.term=[endtimediscountfactor.value '*(' salvagevalue.value ')'];


    case 'discountedsalvagevalue'
        salvagevalue=retrievemodelinformation(ocStruct,'salvagevalue');
        endtimediscountfactor=retrievemodelinformation(ocStruct,'endtimediscountfactor');
        if strcmp(salvagevalue.value,'0')
            algebraic.term=salvagevalue.value;
        else
            algebraic.term=[endtimediscountfactor.value '*(' salvagevalue.value ')'];
        end
    case 'controlconstraint'
        switch classtype
            case 'inequality'
                inequalitycontrolconstraint=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
                algebraic.term=inequalitycontrolconstraint.value;
            otherwise
        end

    case 'constraint'
        inequalitycontrolconstraint=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
        inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
        algebraic.term=[inequalitycontrolconstraint.value inequalitystateconstraint.value];

    case 'symbolicconstraint'
        inequalitycontrolconstraint=retrievemodelinformation(ocStruct,'inequalitycontrolconstraint');
        %inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
        %algebraic.term=[inequalitycontrolconstraint.value inequalitystateconstraint.value];
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numel(inequalitycontrolconstraint.value)
            if ii==1
                algebraic.term{ii}=[symstr '([''[' inequalitycontrolconstraint.value{ii} ';'' ...'];
            elseif ii<numel(inequalitycontrolconstraint.value)
                algebraic.term{ii}=['''' inequalitycontrolconstraint.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' inequalitycontrolconstraint.value{ii} ']''])'];
            end
        end

    case 'stateconstraint'
        switch classtype
            case 'inequality'
                inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
                algebraic.term=inequalitystateconstraint.value;
            otherwise
        end

    case 'canonicalsystem'
        canonicalsystem=retrievemodelinformation(ocStruct,'canonicalsystem','','');
        for ii=1:numel(canonicalsystem.value)
            if ii==1
                algebraic.term{ii}=['[' canonicalsystem.value{ii} '; ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=[canonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=[canonicalsystem.value{ii} ']'];
            end
        end

    case 'canonicalsystemcont'
        canonicalsystem=retrievemodelinformation(ocStruct,'canonicalsystemcont','','');
        for ii=1:numel(canonicalsystem.value)
            if ii==1
                algebraic.term{ii}=['[' canonicalsystem.value{ii} '; ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=[canonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=[canonicalsystem.value{ii} ']'];
            end
        end

    case 'exogenousdynamicsterm'
        exogenousfunctionterm=retrievemodelinformation(ocStruct,'exogenousdynamicsterm');
        exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousdynamicsnum');
        for ii=1:exogenousfunctionnum.value
            if ii==1
                if exogenousfunctionnum.value>1
                    algebraic.term{ii}=['out=[' exogenousfunctionterm.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['out=' exogenousfunctionterm.value{ii}];
                end
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} ']'];
            end
        end

    case 'exogenousfunctionterm'
        exogenousfunctionterm=retrievemodelinformation(ocStruct,'exogenousfunctionterm');
        exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
        for ii=1:exogenousfunctionnum.value
            if ii==1
                if exogenousfunctionnum.value>1
                    algebraic.term{ii}=['out=[' exogenousfunctionterm.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['out=' exogenousfunctionterm.value{ii}];
                end
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} ']'];
            end
        end

    case 'exogenousdynamicsjacobianterm'
        arcidentifier=field2arcidentifier(arcfield);
        exogenousdynamicsjacobianterm=retrievemodelinformation(ocStruct,'exogenousdynamicsjacobianterm',arcidentifier,getsymkernel);
        numterm=numel(exogenousdynamicsjacobianterm.value);
        for ii=1:numterm
            if ii==1 && numterm==1
                algebraic.term{ii}=['[' exogenousdynamicsjacobianterm.value{ii}];
            elseif ii==1 && numterm>1
                algebraic.term{ii}=['[' exogenousdynamicsjacobianterm.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[exogenousdynamicsjacobianterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=exogenousdynamicsjacobianterm.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        
    case 'exogenousdynamicsparameterjacobianterm'
        arcidentifier=field2arcidentifier(arcfield);
        exogenousdynamicsparameterjacobianterm=retrievemodelinformation(ocStruct,'exogenousdynamicsparameterjacobianterm',arcidentifier,getsymkernel);
        numterm=numel(exogenousdynamicsparameterjacobianterm.value);
        for ii=1:numterm
            if ii==1 && numterm==1
                algebraic.term{ii}=['[' exogenousdynamicsparameterjacobianterm.value{ii}];
            elseif ii==1 && numterm>1
                algebraic.term{ii}=['[' exogenousdynamicsparameterjacobianterm.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[exogenousdynamicsparameterjacobianterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=exogenousdynamicsparameterjacobianterm.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'exogenousjacobianterm'
        arcidentifier=field2arcidentifier(arcfield);
        exogenousjacobianterm=retrievemodelinformation(ocStruct,'exogenousjacobianterm',arcidentifier,getsymkernel);
        numterm=numel(exogenousjacobianterm.value);
        for ii=1:numterm
            if ii==1 && numterm==1
                algebraic.term{ii}=['[' exogenousjacobianterm.value{ii}];
            elseif ii==1 && numterm>1
                algebraic.term{ii}=['[' exogenousjacobianterm.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[exogenousjacobianterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=exogenousjacobianterm.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        
    case 'exogenousparameterjacobianterm'
        arcidentifier=field2arcidentifier(arcfield);
        exogenousparameterjacobianterm=retrievemodelinformation(ocStruct,'exogenousparameterjacobianterm',arcidentifier,getsymkernel);
        numterm=numel(exogenousparameterjacobianterm.value);
        for ii=1:numterm
            if ii==1 && numterm==1
                algebraic.term{ii}=['[' exogenousparameterjacobianterm.value{ii}];
            elseif ii==1 && numterm>1
                algebraic.term{ii}=['[' exogenousparameterjacobianterm.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[exogenousparameterjacobianterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=exogenousparameterjacobianterm.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        

    case 'nonsmoothfunctionterm'
        arcidentifier=field2arcidentifier(arcfield);
        nonsmoothfunctionterm=retrievemodelinformation(ocStruct,'nonsmoothfunction',arcidentifier,getsymkernel);
        nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
        for ii=1:nonsmoothfunctionnum.value
            if ii==1
                if nonsmoothfunctionnum.value>1
                    algebraic.term{ii}=['out=[' nonsmoothfunctionterm.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['out=' nonsmoothfunctionterm.value{ii} ';'];
                end
            elseif ii<nonsmoothfunctionnum.value
                algebraic.term{ii}=[nonsmoothfunctionterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=[nonsmoothfunctionterm.value{ii} '];'];
            end
        end

    case 'nonsmoothjacobianterm'
        arcidentifier=field2arcidentifier(arcfield);
        dependentname=getbasicname('dependent');
        nonsmoothjacobianterm=retrievemodelinformation(ocStruct,'nonsmoothjacobianterm',arcidentifier,getsymkernel);
        numterm=numel(nonsmoothjacobianterm.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' 'repmat([' nonsmoothjacobianterm.value{ii} '].'',1,size(' dependentname ',2))'];
            elseif ii<numterm
                algebraic.term{ii}=[nonsmoothjacobianterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=nonsmoothjacobianterm.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'nonsmoothparameterjacobianterm'
        arcidentifier=field2arcidentifier(arcfield);
        nonsmoothparameterjacobianterm=retrievemodelinformation(ocStruct,'nonsmoothparameterjacobianterm',arcidentifier,getsymkernel);
        numterm=numel(nonsmoothparameterjacobianterm.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' nonsmoothparameterjacobianterm.value{ii}];
            elseif ii<numterm
                algebraic.term{ii}=[nonsmoothparameterjacobianterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=nonsmoothparameterjacobianterm.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'switchingvalue'
        switchingvalue=retrievemodelinformation(ocStruct,'switchingvalue',field2arcidentifier(arcfield));
        nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
            for ii=1:nonsmoothfunctionnum.value
                if ii==1
                    if nonsmoothfunctionnum.value>1
                        algebraic.term{ii}=['[' switchingvalue.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' switchingvalue.value{ii}];
                    end
                elseif ii<nonsmoothfunctionnum.value
                    algebraic.term{ii}=[switchingvalue.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[switchingvalue.value{ii}];
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];


    case 'algebraicequation'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequation=retrievemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
        for ii=1:numel(algebraicequation.value)
            if ii==1
                if algebraicequationnum.value>1
                    algebraic.term{ii}=['[' algebraicequation.value{ii} '; ...'];
                else
                    algebraic.term{ii}=algebraicequation.value{ii};
                end
            elseif ii<numel(algebraicequation.value)
                algebraic.term{ii}=['' algebraicequation.value{ii} '; ...'];
            else
                algebraic.term{ii}=['' algebraicequation.value{ii} ']'];
            end
        end

    case 'algebraicequationcont'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequation=retrievemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
        for ii=1:numel(algebraicequation.value)
            if ii==1
                if algebraicequationnum.value>1
                    algebraic.term{ii}=['contpar*[' algebraicequation.value{ii} '; ...'];
                else
                    algebraic.term{ii}=algebraicequation.value{ii};
                end
            elseif ii<numel(algebraicequation.value)
                algebraic.term{ii}=['' algebraicequation.value{ii} '; ...'];
            else
                algebraic.term{ii}=['' algebraicequation.value{ii} ']'];
            end
        end

    case 'algebraicequationimplicit'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequationimplicit=retrievemodelinformation(ocStruct,'algebraicequationimplicit',arcidentifier,getsymkernel);
        algebraicequationimplicitnum=numel(algebraicequationimplicit.value);
        for ii=1:algebraicequationimplicitnum
            if ii==1
                if algebraicequationimplicitnum>1
                    algebraic.term{ii}=['[' algebraicequationimplicit.value{ii} '; ...'];
                else
                    if isempty(algebraicequationimplicit.value{ii})
                        algebraic.term{ii}='[]';
                    else
                        algebraic.term{ii}=algebraicequationimplicit.value{ii};
                    end
                end
            elseif ii<numel(algebraicequationimplicit.value)
                algebraic.term{ii}=['' algebraicequationimplicit.value{ii} '; ...'];
            else
                algebraic.term{ii}=['' algebraicequationimplicit.value{ii} ']'];
            end
        end

    case 'algebraicequationimplicitjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequationimplicitjacobian=retrievemodelinformation(ocStruct,'algebraicequationimplicitjacobian',arcidentifier,getsymkernel);
        numterm=numel(algebraicequationimplicitjacobian.value);
        for ii=1:numterm
            if ii==1 && numterm==1
                if strcmp(algebraicequationimplicitjacobian.value,'[]')
                    algebraic.term{ii}=[algebraicequationimplicitjacobian.value{ii}];
                else
                    algebraic.term{ii}=['[' algebraicequationimplicitjacobian.value{ii} ']'];
                end
            elseif ii==1 && numterm>1
                algebraic.term{ii}=['[' algebraicequationimplicitjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[algebraicequationimplicitjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=algebraicequationimplicitjacobian.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'algebraicequationimplicitparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequationimplicitparameterjacobian=retrievemodelinformation(ocStruct,'algebraicequationimplicitparameterjacobian',arcidentifier,getsymkernel);
        numterm=numel(algebraicequationimplicitparameterjacobian.value);
        for ii=1:numterm
            if ii==1 && numterm==1
                if strcmp(algebraicequationimplicitparameterjacobian.value,'[]')
                    algebraic.term{ii}=[algebraicequationimplicitparameterjacobian.value{ii}];
                else
                    algebraic.term{ii}=['[' algebraicequationimplicitparameterjacobian.value{ii} ']'];
                end
            elseif ii==1 && numterm>1
                algebraic.term{ii}=['[' algebraicequationimplicitparameterjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[algebraicequationimplicitparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=algebraicequationimplicitparameterjacobian.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'symboliccanonicalsystem'
        canonicalsystem=retrievemodelinformation(ocStruct,'canonicalsystem','','');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numel(canonicalsystem.value)
            if ii==1
                algebraic.term{ii}=[symstr '([''[' canonicalsystem.value{ii} ';'' ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=['''' canonicalsystem.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystem.value{ii} ']''])'];
            end
        end

    case 'symbolicalgebraicequation'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequation=retrievemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numel(algebraicequation.value)
            if ii==1
                if algebraicequationnum.value>1
                    algebraic.term{ii}=[symstr '([''[' algebraicequation.value{ii} ';'' ...'];
                else
                    algebraic.term{ii}=[symstr '(''' algebraicequation.value{ii} ''')'];
                end
            elseif ii<numel(algebraicequation.value)
                algebraic.term{ii}=['''' algebraicequation.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' algebraicequation.value{ii} ']''])'];
            end
        end

    case 'symboliccanonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrievemodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=[symstr '([''[' canonicalsystemjacobian.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ']'''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'symbolicstatecontrolcanonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrievemodelinformation(ocStruct,'canonicalsystemjacobianstatecontrol',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=[symstr '([''' canonicalsystemjacobian.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' canonicalsystemjacobian.value{ii} ''''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'symbolicstatecontrolcanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsysteminstatecontrol=retrievemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,getsymkernel);
        numterm=numel(specificcanonicalsysteminstatecontrol.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=[symstr '([''[' specificcanonicalsysteminstatecontrol.value{ii} ';'' ...'];
            elseif ii<numterm
                algebraic.term{ii}=['''' specificcanonicalsysteminstatecontrol.value{ii} ';'' ...'];
            else
                algebraic.term{ii}=['''' specificcanonicalsysteminstatecontrol.value{ii} ']'''];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} '])'];

    case 'specificcanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsystem=retrievemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,getsymkernel);
        numterm=numel(specificcanonicalsystem.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' specificcanonicalsystem.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[specificcanonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=specificcanonicalsystem.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'specificpontryaginfunction'
        arcidentifier=field2arcidentifier(arcfield);
        specificpontryaginfunction=retrievemodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,getsymkernel);
        
        algebraic.term=specificpontryaginfunction.value;

    case 'symboliccanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsystem=retrievemodelinformation(ocStruct,'symboliccanonicalsystem',arcidentifier,getsymkernel);
        numterm=numel(specificcanonicalsystem.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' specificcanonicalsystem.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[specificcanonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=specificcanonicalsystem.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'statecostatejacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrievemodelinformation(ocStruct,'statecostatejacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobian=retrievemodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemjacobianstatecontrol'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobianstatecontrol=retrievemodelinformation(ocStruct,'canonicalsystemjacobianstatecontrol',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobianstatecontrol.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobianstatecontrol.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobianstatecontrol.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobianstatecontrol.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'optimalcontroldynamicsleftside'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'optimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsleftside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end

    case 'optimalcontroldynamics'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'optimalcontroldynamics',arcidentifier,getsymkernel);
        %implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex',arcidentifier);
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsleftside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                end
            end
            if numterm>1
                algebraic.term{ii}=[algebraic.term{ii} ']'];
            end
        else
            algebraic.term=[];
        end

    case 'symbolicoptimalcontroldynamics'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'optimalcontroldynamics',arcidentifier,getsymkernel);
        %implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex',arcidentifier);
        numterm=numel(optimalcontroldynamicsleftside.value);
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=[symstr '([''[' optimalcontroldynamicsleftside.value{ii} ';'' ...'];
                    else
                        algebraic.term{ii}=[symstr '(''' optimalcontroldynamicsleftside.value{ii} ''')'];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=['''' optimalcontroldynamicsleftside.value{ii} ';'' ...'];
                else
                    algebraic.term{ii}=['''' optimalcontroldynamicsleftside.value{ii} ']''])'];
                end
            end
            if numterm>1
                algebraic.term{ii}=[algebraic.term{ii} ']'];
            end
        else
            algebraic.term=[];
        end

    case 'invoptimalcontroldynamicsleftside'
        arcidentifier=field2arcidentifier(arcfield);
        optimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        numterm=numel(optimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' optimalcontroldynamicsleftside.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsleftside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsleftside.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end

    case 'tensoroptimalcontroldynamicsleftside'
        arcidentifier=field2arcidentifier(arcfield);
        tensoroptimalcontroldynamicsleftside=retrievemodelinformation(ocStruct,'tensoroptimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        nummat=numel(tensoroptimalcontroldynamicsleftside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            counter=0;
            for jj=1:nummat
                numterm=numel(tensoroptimalcontroldynamicsleftside.value{jj});
                for ii=1:numterm
                    counter=counter+1;
                    if ii==1
                        if numterm>1
                            algebraic.term{counter}=['out(:,:,' num2str(jj) ')=[' tensoroptimalcontroldynamicsleftside.value{jj}{ii} '; ...'];
                        else
                            algebraic.term{counter}=['out(:,:,' num2str(jj) ')=[' tensoroptimalcontroldynamicsleftside.value{jj}{ii}];
                        end
                    elseif ii<numterm
                        algebraic.term{counter}=[tensoroptimalcontroldynamicsleftside.value{jj}{ii} '; ...'];
                    else
                        algebraic.term{counter}=tensoroptimalcontroldynamicsleftside.value{jj}{ii};
                    end
                end
                algebraic.term{counter}=[algebraic.term{counter} '];'];
            end
        else
            algebraic.term=[];
        end


    case 'jacobianoptimalcontroldynamicsrightside'
        arcidentifier=field2arcidentifier(arcfield);
        jacobianoptimalcontroldynamicsrightside=retrievemodelinformation(ocStruct,'jacobianoptimalcontroldynamicsrightside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        numterm=numel(jacobianoptimalcontroldynamicsrightside.value);
        if ~isempty(implicitnonlinearcontrolindex.value)
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' jacobianoptimalcontroldynamicsrightside.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' jacobianoptimalcontroldynamicsrightside.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[jacobianoptimalcontroldynamicsrightside.value{ii} '; ...'];
                else
                    algebraic.term{ii}=jacobianoptimalcontroldynamicsrightside.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end

    case 'pontryaginfunctionDuDX'
        pontryaginfunctionDuDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDX','','');
        for ii=1:numel(pontryaginfunctionDuDX.value)
            if ii==1
                if numel(pontryaginfunctionDuDX.value)>1
                    algebraic.term{ii}=['[' pontryaginfunctionDuDX.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['[' pontryaginfunctionDuDX.value{ii}];
                end
            elseif ii<numel(pontryaginfunctionDuDX.value)
                algebraic.term{ii}=[pontryaginfunctionDuDX.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDuDX.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDlmmcDX'
        arcidentifier=field2arcidentifier(arcfield);
        pontryaginfunctionDlmmcDX=retrievemodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,getsymkernel);
        if ~isempty(pontryaginfunctionDlmmcDX.value)
            for ii=1:numel(pontryaginfunctionDlmmcDX.value)
                if ii==1
                    if numel(pontryaginfunctionDlmmcDX.value)>1
                        algebraic.term{ii}=['[' pontryaginfunctionDlmmcDX.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=['[' pontryaginfunctionDlmmcDX.value{ii}];
                    end
                elseif ii<numel(pontryaginfunctionDlmmcDX.value)
                    algebraic.term{ii}=[pontryaginfunctionDlmmcDX.value{ii} '; ...'];
                else
                    algebraic.term{ii}=pontryaginfunctionDlmmcDX.value{ii};
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        else
            algebraic.term=[];
        end


    case 'canonicalsystemparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemparameterjacobian=retrievemodelinformation(ocStruct,'canonicalsystemparameterjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemparameterjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemparameterjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];


    case 'statecontrolcanonicalsystemparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        statecontrolcanonicalsystemparameterjacobian=retrievemodelinformation(ocStruct,'statecontrolcanonicalsystemparameterjacobian',arcidentifier,getsymkernel);
        numterm=numel(statecontrolcanonicalsystemparameterjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' statecontrolcanonicalsystemparameterjacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[statecontrolcanonicalsystemparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=statecontrolcanonicalsystemparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemderivativetime'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemderivativetime=retrievemodelinformation(ocStruct,'canonicalsystemderivativetime',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemderivativetime.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemderivativetime.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemderivativetime.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemderivativetime.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemhessian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemhessian=retrievemodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier,getsymkernel);
        nummatrix=numel(canonicalsystemhessian.value);
        numlines=numel(canonicalsystemhessian.value{1});
        counter=0;
        for ii=1:nummatrix
            for jj=1:numlines
                counter=counter+1;
                if jj==1
                    algebraic.term{counter}=['[' canonicalsystemhessian.value{ii}{jj} '; ...'];
                elseif jj<numlines
                    algebraic.term{counter}=[canonicalsystemhessian.value{ii}{jj} '; ...'];
                elseif jj==numlines
                    algebraic.term{counter}=canonicalsystemhessian.value{ii}{jj};
                end
            end
            algebraic.term{counter}=[algebraic.term{counter} ']'];
        end
        algebraic.info=numlines; % return number of matrices for Hessian

    case 'canonicalsystemtotalhessian'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemparameterhessian=retrievemodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier,getsymkernel);
        nummatrix=numel(canonicalsystemparameterhessian.value);
        numlines=numel(canonicalsystemparameterhessian.value{1});
        counter=0;
        for ii=1:nummatrix
            for jj=1:numlines
                counter=counter+1;
                if jj==1
                    algebraic.term{counter}=['[' canonicalsystemparameterhessian.value{ii}{jj} '; ...'];
                elseif jj<numlines
                    algebraic.term{counter}=[canonicalsystemparameterhessian.value{ii}{jj} '; ...'];
                elseif jj==numlines
                    algebraic.term{counter}=canonicalsystemparameterhessian.value{ii}{jj};
                end
            end
            algebraic.term{counter}=[algebraic.term{counter} ']'];
        end
        algebraic.info=numlines; % return number of matrices for Hessian


    case 'canonicalsystemvariationparameterderivative'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemvariationparameterderivative=retrievemodelinformation(ocStruct,'canonicalsystemvariationparameterderivative',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemvariationparameterderivative.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemvariationparameterderivative.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemvariationparameterderivative.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemvariationparameterderivative.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        
    case 'parametervariables'
        parametername=retrievemodelinformation(ocStruct,'parametername');
        algebraic.term=parametername.value;

    case 'discobjectivefunction'
        discobjectivefunction=retrievemodelinformation(ocStruct,'discobjectivefunction');
        algebraic.term{1}=discobjectivefunction.value;

    case 'discobjectivefunctionderivativetime'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionderivativetime=retrievemodelinformation(ocStruct,'discobjectivefunctionderivativetime',arcidentifier,getsymkernel);
        algebraic.term=discobjectivefunctionderivativetime.value;


    case 'symbolicdiscobjectivefunction'
        discobjectivefunction=retrievemodelinformation(ocStruct,'discobjectivefunction');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        algebraic.term{1}=[symstr '(''' discobjectivefunction.value ''')'];

    case 'discobjectivefunctionjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrievemodelinformation(ocStruct,'discobjectivefunctionDX',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' discobjectivefunctionjacobian.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii} ', ...'];
            else
                algebraic.term{ii}=discobjectivefunctionjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'discobjectivefunctionparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrievemodelinformation(ocStruct,'discobjectivefunctionDP',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' discobjectivefunctionjacobian.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[discobjectivefunctionjacobian.value{ii} ', ...'];
            else
                algebraic.term{ii}=discobjectivefunctionjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'hamiltonianfunction'
        hamiltonianfunction=retrievemodelinformation(ocStruct,'hamiltonianfunction');
        algebraic.term{1}=hamiltonianfunction.value;

    case 'pontryaginfunction'
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
        algebraic.term{1}=pontryaginfunction.value;

    case 'symbolicpontryaginfunction'
        pontryaginfunction=retrievemodelinformation(ocStruct,'pontryaginfunction');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        algebraic.term{1}=[symstr '(''' pontryaginfunction.value ''')'];

    case 'pontryaginfunctionDx'
        pontryaginfunctionDx=retrievemodelinformation(ocStruct,'pontryaginfunctionDx');
        numterm=numel(pontryaginfunctionDx.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDx.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDx.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDx.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'implicitcontroldynamicsGX'
        implicitcontroldynamicsGX=retrievemodelinformation(ocStruct,'implicitcontroldynamicsGX',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(implicitcontroldynamicsGX.value);
        for ii=1:numterm
            if ii==1
                if numterm==1
                    algebraic.term{ii}=['[' implicitcontroldynamicsGX.value{ii}];
                else
                    algebraic.term{ii}=['[' implicitcontroldynamicsGX.value{ii} '; ...'];
                end
            elseif ii<numterm
                algebraic.term{ii}=[implicitcontroldynamicsGX.value{ii} '; ...'];
            else
                algebraic.term{ii}=implicitcontroldynamicsGX.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'implicitcontroldynamicsfunction'
        implicitcontroldynamicsfunction=retrievemodelinformation(ocStruct,'implicitcontroldynamicsfunction',field2arcidentifier(arcfield),getsymkernel);
        if ~isempty(implicitcontroldynamicsfunction.value)
        symimplicitcontroldynamicsfunction=mystr2sym(implicitcontroldynamicsfunction.value);
        implicitcontroldynamicsfunction.value=string2cell(removematrixstring(char(symimplicitcontroldynamicsfunction(:))),'matrix');
        else
            implicitcontroldynamicsfunction.value{1}='[]';
        end
        numterm=numel(implicitcontroldynamicsfunction.value);
        for ii=1:numterm
            if ii==1
                if numterm==1
                    algebraic.term{ii}=[implicitcontroldynamicsfunction.value{ii}];
                else
                    algebraic.term{ii}=['[' implicitcontroldynamicsfunction.value{ii} '; ...'];
                end
            elseif ii<numterm
                algebraic.term{ii}=[implicitcontroldynamicsfunction.value{ii} '; ...'];
            else
                algebraic.term{ii}=implicitcontroldynamicsfunction.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'implicitcontroldynamicsGP'
        implicitcontroldynamicsGP=retrievemodelinformation(ocStruct,'implicitcontroldynamicsGP',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(implicitcontroldynamicsGP.value);
        for ii=1:numterm
            if ii==1
                if numterm==1
                    algebraic.term{ii}=['[' implicitcontroldynamicsGP.value{ii}];
                else
                    algebraic.term{ii}=['[' implicitcontroldynamicsGP.value{ii} '; ...'];
                end
            elseif ii<numterm
                algebraic.term{ii}=[implicitcontroldynamicsGP.value{ii} '; ...'];
            else
                algebraic.term{ii}=implicitcontroldynamicsGP.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'implicitcontroldynamicsGU'
        implicitcontroldynamicsGU=retrievemodelinformation(ocStruct,'implicitcontroldynamicsGU',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(implicitcontroldynamicsGU.value);
        for ii=1:numterm
            if ii==1
                if numterm==1
                    algebraic.term{ii}=['[' implicitcontroldynamicsGU.value{ii}];
                else
                    algebraic.term{ii}=['[' implicitcontroldynamicsGU.value{ii} '; ...'];
                end
            elseif ii<numterm
                algebraic.term{ii}=[implicitcontroldynamicsGU.value{ii} '; ...'];
            else
                algebraic.term{ii}=implicitcontroldynamicsGU.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'optimalcontroldynamicsPuu'
        optimalcontroldynamicsPuu=retrievemodelinformation(ocStruct,'optimalcontroldynamicsPuu',field2arcidentifier(arcfield),getsymkernel);
        %totalvariable=[controlname.value lagrangemultipliercontrolname.value];
        numterm=numel(optimalcontroldynamicsPuu.value);
        for ii=1:numterm
            if numterm==1
                %if length(totalvariable)==1
                    algebraic.term{ii}=optimalcontroldynamicsPuu.value{ii};
                %else
                %    algebraic.term{ii}=['[' pontryaginfunctionDU2.value{ii}];
                %end
            else
                if ii==1
                    algebraic.term{ii}=['[' optimalcontroldynamicsPuu.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[optimalcontroldynamicsPuu.value{ii} '; ...'];
                else
                    algebraic.term{ii}=optimalcontroldynamicsPuu.value{ii};
                end
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
        
    case 'dlagrangiandU'
        implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',field2arcidentifier(arcfield));
        dlagrangiandU=retrievemodelinformation(ocStruct,'dlagrangiandU',field2arcidentifier(arcfield),getsymkernel);
        if ~isempty(implicitnonlinearcontrolindex.value)
            dlagrangiandU.value=dlagrangiandU.value(implicitnonlinearcontrolindex.value);
            dlagrangiandU.value=string2cell(removematrixstring(char(dlagrangiandU.value)),'matrix');
            zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',field2arcidentifier(arcfield));
            lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
            lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        else
            dlagrangiandU.value=[];
            dlagrangiandU.value{1}='[]';
        end
        numterm=numel(dlagrangiandU.value);
        for ii=1:numterm
            if numterm==1
                algebraic.term{ii}=dlagrangiandU.value{ii};
            else
                if ii==1
                    algebraic.term{ii}=['[' dlagrangiandU.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[dlagrangiandU.value{ii} '; ...'];
                else
                    algebraic.term{ii}=dlagrangiandU.value{ii};
                end
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'd2lagrangiandU2'
        d2lagrangiandU2=retrievemodelinformation(ocStruct,'d2lagrangiandU2',field2arcidentifier(arcfield),getsymkernel);
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',field2arcidentifier(arcfield));
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        %totalvariable=[controlname.value lagrangemultipliercontrolname.value];
        numterm=numel(d2lagrangiandU2.value);
        for ii=1:numterm
            if numterm==1
                %if length(totalvariable)==1
                    algebraic.term{ii}=d2lagrangiandU2.value{ii};
                %else
                %    algebraic.term{ii}=['[' d2lagrangiandU2.value{ii}];
                %end
            else
                if ii==1
                    algebraic.term{ii}=['[' d2lagrangiandU2.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[d2lagrangiandU2.value{ii} '; ...'];
                else
                    algebraic.term{ii}=d2lagrangiandU2.value{ii};
                end
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'd2lagrangiandui2'
        d2lagrangiandui2=retrievemodelinformation(ocStruct,'d2lagrangiandui2',field2arcidentifier(arcfield),getsymkernel);
        zerolmmcindex=retrievemodelinformation(ocStruct,'zerolmmcindex',field2arcidentifier(arcfield));
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        lagrangemultipliercontrolname.value(zerolmmcindex.value)=[];
        %totalvariable=[controlname.value lagrangemultipliercontrolname.value];
        numterm=numel(d2lagrangiandui2.value);
        for ii=1:numterm
            if numterm==1
                %if length(totalvariable)==1
                    algebraic.term{ii}=d2lagrangiandui2.value{ii};
                %else
                %    algebraic.term{ii}=['[' d2lagrangiandui2.value{ii}];
                %end
            else
                if ii==1
                    algebraic.term{ii}=['[' d2lagrangiandui2.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[d2lagrangiandui2.value{ii} '; ...'];
                else
                    algebraic.term{ii}=d2lagrangiandui2.value{ii};
                end
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'd2lagrangiandUdX'
        d2lagrangiandUdX=retrievemodelinformation(ocStruct,'d2lagrangiandUdX',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(d2lagrangiandUdX.value);
        for ii=1:numterm
            if numterm==1
                if strcmp(d2lagrangiandUdX.value{ii},'[]')
                    algebraic.term{ii}='[';
                else
                    algebraic.term{ii}=['[' d2lagrangiandUdX.value{ii}];
                end
            else
                if ii==1
                    algebraic.term{ii}=['[' d2lagrangiandUdX.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[d2lagrangiandUdX.value{ii} '; ...'];
                else
                    algebraic.term{ii}=d2lagrangiandUdX.value{ii};
                end
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'D2pontryaginfunctionDUDX'
        D2pontryaginfunctionDUDX=retrievemodelinformation(ocStruct,'D2pontryaginfunctionDUDX',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(D2pontryaginfunctionDUDX.value);
        for ii=1:numterm
            if numterm==1
                if strcmp(D2pontryaginfunctionDUDX.value{ii},'[]')
                    algebraic.term{ii}='[';
                else
                    algebraic.term{ii}=['[' D2pontryaginfunctionDUDX.value{ii}];
                end
            else
                if ii==1
                    algebraic.term{ii}=['[' D2pontryaginfunctionDUDX.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[D2pontryaginfunctionDUDX.value{ii} '; ...'];
                else
                    algebraic.term{ii}=D2pontryaginfunctionDUDX.value{ii};
                end
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'D2pontryaginfunctionDui2'
        D2pontryaginfunctionDui2=retrievemodelinformation(ocStruct,'D2pontryaginfunctionDui2',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(D2pontryaginfunctionDui2.value);
        for ii=1:numterm
            if numterm==1
                if strcmp(D2pontryaginfunctionDui2.value{ii},'[]')
                    algebraic.term{ii}='[';
                else
                    algebraic.term{ii}=['[' D2pontryaginfunctionDui2.value{ii}];
                end
            else
                if ii==1
                    algebraic.term{ii}=['[' D2pontryaginfunctionDui2.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[D2pontryaginfunctionDui2.value{ii} '; ...'];
                else
                    algebraic.term{ii}=D2pontryaginfunctionDui2.value{ii};
                end
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'd2lagrangiandUdP'
        d2lagrangiandUdP=retrievemodelinformation(ocStruct,'d2lagrangiandUdP',field2arcidentifier(arcfield),getsymkernel);
        numterm=numel(d2lagrangiandUdP.value);
        for ii=1:numterm
            if numterm==1
                if strcmp(d2lagrangiandUdP.value{ii},'[]')
                    algebraic.term{ii}='[';
                else
                    algebraic.term{ii}=['[' d2lagrangiandUdP.value{ii}];
                end
            else
                if ii==1
                    algebraic.term{ii}=['[' d2lagrangiandUdP.value{ii} '; ...'];
                elseif ii<numterm
                    algebraic.term{ii}=[d2lagrangiandUdP.value{ii} '; ...'];
                else
                    algebraic.term{ii}=d2lagrangiandUdP.value{ii};
                end
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDu2'
        pontryaginfunctionDu2=retrievemodelinformation(ocStruct,'pontryaginfunctionDu2',[],getsymkernel);
        numterm=numel(pontryaginfunctionDu2.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDu2.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDu2.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDu2.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];


    case 'D2pontryaginfunctionDX2'
        D2pontryaginfunctionDX2=retrievemodelinformation(ocStruct,'D2pontryaginfunctionDX2',[],getsymkernel);
        numterm=numel(D2pontryaginfunctionDX2.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' D2pontryaginfunctionDX2.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[D2pontryaginfunctionDX2.value{ii} '; ...'];
            else
                algebraic.term{ii}=D2pontryaginfunctionDX2.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDuDx'
        pontryaginfunctionDuDx=retrievemodelinformation(ocStruct,'pontryaginfunctionDuDx');
        numterm=numel(pontryaginfunctionDuDx.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDuDx.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDuDx.value{ii} '; ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDuDx.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'pontryaginfunctionDX'
        arcidentifier=field2arcidentifier(arcfield);
        pontryaginfunctionDx=retrievemodelinformation(ocStruct,'pontryaginfunctionDX',arcidentifier,getsymkernel);
        numterm=numel(pontryaginfunctionDx.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' pontryaginfunctionDx.value{ii} ', ...'];
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDx.value{ii} ', ...'];
            else
                algebraic.term{ii}=pontryaginfunctionDx.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'DpontryaginfunctionDu'
        arcidentifier=field2arcidentifier(arcfield);
        DpontryaginfunctionDu=retrievemodelinformation(ocStruct,'DpontryaginfunctionDu',arcidentifier,getsymkernel);
        if ~isempty(DpontryaginfunctionDu.value{1})
            numterm=numel(DpontryaginfunctionDu.value);
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' DpontryaginfunctionDu.value{ii} ', ...'];
                    else
                        algebraic.term{ii}=[DpontryaginfunctionDu.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[DpontryaginfunctionDu.value{ii} ', ...'];
                else
                    algebraic.term{ii}=DpontryaginfunctionDu.value{ii};
                end
            end
            if numterm>1
                algebraic.term{ii}=[algebraic.term{ii} ']'];
            end
        else
            algebraic.term{1}='[]';
        end

    case 'DHamiltonianDu'
        arcidentifier=field2arcidentifier(arcfield);
        DHamiltonianDu=retrievemodelinformation(ocStruct,'DHamiltonianDu',arcidentifier,getsymkernel);
        if ~isempty(DHamiltonianDu.value{1})
            numterm=numel(DHamiltonianDu.value);
            for ii=1:numterm
                if ii==1
                    if numterm>1
                        algebraic.term{ii}=['[' DHamiltonianDu.value{ii} '; ...'];
                    else
                        algebraic.term{ii}=[DHamiltonianDu.value{ii}];
                    end
                elseif ii<numterm
                    algebraic.term{ii}=[DHamiltonianDu.value{ii} '; ...'];
                else
                    algebraic.term{ii}=DHamiltonianDu.value{ii};
                end
            end
            if numterm>1
                algebraic.term{ii}=[algebraic.term{ii} ']'];
            end
        else
            algebraic.term{1}='[]';
        end

    case 'inequalitystateconstrainttimederivative'
        inequalitystateconstraint=retrievemodelinformation(ocStruct,'inequalitystateconstraint');
        inequalitystateconstrainttimederivative=retrievemodelinformation(ocStruct,'inequalitystateconstrainttimederivative');
        inequalitystateconstraintorder=retrievemodelinformation(ocStruct,'inequalitystateconstraintorder');
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        algebraic.term{1}='';
        for jj=1:inequalitystateconstraintnum.value
            algebraic.term{1}=[algebraic.term{1} inequalitystateconstraint.value{jj} ';'];
        end
        algebraic.term{1}(end)=[];
        for ii=(1:inequalitystateconstraintorder.value)+1
            algebraic.term{ii}='';
            for jj=1:inequalitystateconstraintnum.value
                algebraic.term{ii}=[algebraic.term{ii} inequalitystateconstrainttimederivative.value{jj}{ii-1} ';'];
            end
            algebraic.term{ii}(end)='';
        end


    case 'hamiltonianDu'
        hamiltonianDu=retrievemodelinformation(ocStruct,'hamiltonianDu','',getsymkernel);
        numterm=numel(hamiltonianDu.value);
        for ii=1:numterm
            if ii==1
                if numterm>1
                    algebraic.term{ii}=['[' hamiltonianDu.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[hamiltonianDu.value{ii}];
                end
            elseif ii<numterm
                algebraic.term{ii}=[hamiltonianDu.value{ii} '; ...'];
            else
                algebraic.term{ii}=hamiltonianDu.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end

    case 'statedynamics'
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics','','');
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[statedynamics.value{ii}];
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end

    case 'statedynamics'
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics','','');
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[statedynamics.value{ii}];
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end
        
    case 'adjointsystem'
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','0','');
        for ii=1:numel(adjointsystem.value)
            if ii==1
                if numel(adjointsystem.value)>1
                    algebraic.term{ii}=['[' adjointsystem.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[adjointsystem.value{ii}];
                end
            elseif ii<numel(adjointsystem.value)
                algebraic.term{ii}=[adjointsystem.value{ii} '; ...'];
            else
                if numel(adjointsystem.value)>1
                    algebraic.term{ii}=[adjointsystem.value{ii} ']'];
                end
            end
        end

    case 'objectivefunction'
        objectivefunction=retrievemodelinformation(ocStruct,'objectiveintegrand');
        algebraic.term{1}=objectivefunction.value;
        
    case 'transversalitycondition'
        transversalitycondition=retrievemodelinformation(ocStruct,'transversalitycondition','',getsymkernel);
        for ii=1:numel(transversalitycondition.value)
            if ii==1
                if numel(transversalitycondition.value)>1
                    algebraic.term{ii}=['[' transversalitycondition.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[transversalitycondition.value{ii}];
                end
            elseif ii<numel(transversalitycondition.value)
                algebraic.term{ii}=[transversalitycondition.value{ii} '; ...'];
            else
                if numel(transversalitycondition.value)>1
                    algebraic.term{ii}=[transversalitycondition.value{ii} ']'];
                end
            end
        end
        
    case 'constraintjacobian'
        constraintjacobian=retrievemodelinformation(ocStruct,'constraintjacobian','',getsymkernel);
        for ii=1:numel(constraintjacobian.value)
            if ii==1
                if numel(constraintjacobian.value)>1
                    algebraic.term{ii}=['[' constraintjacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[constraintjacobian.value{ii}];
                end
            elseif ii<numel(constraintjacobian.value)
                algebraic.term{ii}=[constraintjacobian.value{ii} '; ...'];
            else
                if numel(constraintjacobian.value)>1
                    algebraic.term{ii}=[constraintjacobian.value{ii} '];'];
                end
            end
        end
        
    case 'constraintparameterjacobian'
        constraintjacobian=retrievemodelinformation(ocStruct,'constraintparameterjacobian','',getsymkernel);
        for ii=1:numel(constraintjacobian.value)
            if ii==1
                if numel(constraintjacobian.value)>1
                    algebraic.term{ii}=['[' constraintjacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[constraintjacobian.value{ii}];
                end
            elseif ii<numel(constraintjacobian.value)
                algebraic.term{ii}=[constraintjacobian.value{ii} '; ...'];
            else
                if numel(constraintjacobian.value)>1
                    algebraic.term{ii}=[constraintjacobian.value{ii} '];'];
                end
            end
        end
        
    case 'constraintcontroljacobian'
        constraintjacobian=retrievemodelinformation(ocStruct,'constraintcontroljacobian','',getsymkernel);
        for ii=1:numel(constraintjacobian.value)
            if ii==1
                if numel(constraintjacobian.value)>1
                    algebraic.term{ii}=['[' constraintjacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[constraintjacobian.value{ii}];
                end
            elseif ii<numel(constraintjacobian.value)
                algebraic.term{ii}=[constraintjacobian.value{ii} '; ...'];
            else
                if numel(constraintjacobian.value)>1
                    algebraic.term{ii}=[constraintjacobian.value{ii} '];'];
                end
            end
        end
        
    case 'variationalcontrolvalue'
        variationaloptimalvalue=retrievemodelinformation(ocStruct,'variationaloptimalvalue',field2arcidentifier(arcfield),getsymkernel);
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        if ~controlnum.value
            algebraic.term{1}='[]';
        else
            for ii=1:controlnum.value
                if ii==1
                    if controlnum.value>1
                        algebraic.term{ii}=['[' variationaloptimalvalue.value.(controlname.value{ii}) '; ...'];
                    else
                        algebraic.term{ii}=['[' variationaloptimalvalue.value.(controlname.value{ii})];
                    end
                elseif ii<controlnum.value
                    algebraic.term{ii}=[variationaloptimalvalue.value.(controlname.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=[variationaloptimalvalue.value.(controlname.value{ii})];
                end
            end
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
        
    case 'admissiblefactor'
        constraintcombinationindex=retrievemodelinformation(ocStruct,'constraintcombinationindex',field2arcidentifier(arcfield));
        cstridx=find(~constraintcombinationindex.value);
        lmcidx=find(constraintcombinationindex.value);
        lagmcc=getbasicname('lagrangemultcc');
        cstr=getbasicname('constraint');
        if isempty(cstridx)
            cstridxchar='';
        else
            cstridxchar=strrep(num2str(cstridx),'  ',' ');
            if length(cstridx)>1
                cstridxchar=['[' cstridxchar ']'];
            end
        end
        if isempty(lmcidx)
            lmcidxchar='';
        else
            lmcidxchar=strrep(num2str(lmcidx),'  ',' ');
            if length(lmcidxchar)>1
                lmcidxchar=['[' lmcidxchar ']'];
            end
        end
        if isempty(cstridxchar) && isempty(lmcidxchar)
            algebraic.term='[]';
        elseif isempty(cstridxchar) && ~isempty(lmcidxchar)
            algebraic.term=sprintf('prod(%s(%s))/maxfac',lagmcc,lmcidxchar);
        elseif ~isempty(cstridxchar) && isempty(lmcidxchar)
            algebraic.term=sprintf('min(%s(%s))',cstr,cstridxchar);
        else
            algebraic.term=sprintf('min([%s(%s);%s(%s)])',cstr, cstridxchar,lagmcc,lmcidxchar);
        end
        data.type='char';
        data.description='derivative of the implicitly give generalized control dynamics with respect to the generalized state.';
        
    case 'variationalparameterobjectivefunction'
        variationalparameterobjectivefunction=retrievemodelinformation(ocStruct,'variationalparameterobjectivefunction',field2arcidentifier(arcfield),getsymkernel);
        algebraic.term{1}=variationalparameterobjectivefunction.value;
    case  'variationaljacobian'
        arcidentifier=field2arcidentifier(arcfield);
        variationaljacobian=retrievemodelinformation(ocStruct,'variationaljacobian',arcidentifier,getsymkernel);
        numterm=numel(variationaljacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' variationaljacobian.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[variationaljacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=variationaljacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];
        
    case  'variationalhamiltonianfunction'
        arcidentifier=field2arcidentifier(arcfield);
        variationalhamiltonianfunction=retrievemodelinformation(ocStruct,'variationalhamiltonianfunction',arcidentifier,getsymkernel);
        algebraic.term{1}=variationalhamiltonianfunction.value;

    case  'parameterderivativehamiltonianfunction'
        arcidentifier=field2arcidentifier(arcfield);
        parameterderivativehamiltonianfunction=retrievemodelinformation(ocStruct,'parameterderivativehamiltonianfunction',arcidentifier,getsymkernel);
        algebraic.term{1}=parameterderivativehamiltonianfunction.value;

        
    case  'variationalsalvagevalue'
        variationalsalvagevalue=retrievemodelinformation(ocStruct,'variationalsalvagevalue',[],getsymkernel);
        for ii=1:numel(variationalsalvagevalue.value)
            if ii==1
                if numel(variationalsalvagevalue.value)>1
                    algebraic.term{ii}=['[' variationalsalvagevalue.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[variationalsalvagevalue.value{ii}];
                end
            elseif ii<numel(variationalsalvagevalue.value)
                algebraic.term{ii}=[variationalsalvagevalue.value{ii} '; ...'];
            else
                if numel(variationalsalvagevalue.value)>1
                    algebraic.term{ii}=[variationalsalvagevalue.value{ii} ']'];
                end
            end
        end
        
    case  'parametervariationalsalvagevalue'
        parametervariationalsalvagevalue=retrievemodelinformation(ocStruct,'parametervariationalsalvagevalue',[],getsymkernel);
        for ii=1:numel(parametervariationalsalvagevalue.value)
            if ii==1
                if numel(parametervariationalsalvagevalue.value)>1
                    algebraic.term{ii}=['[' parametervariationalsalvagevalue.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[parametervariationalsalvagevalue.value{ii}];
                end
            elseif ii<numel(parametervariationalsalvagevalue.value)
                algebraic.term{ii}=[parametervariationalsalvagevalue.value{ii} '; ...'];
            else
                if numel(parametervariationalsalvagevalue.value)>1
                    algebraic.term{ii}=[parametervariationalsalvagevalue.value{ii} ']'];
                end
            end
        end
        
    case 'variationaltransversalitycondition'
        variationaltransversalitycondition=retrievemodelinformation(ocStruct,'variationaltransversalitycondition','',getsymkernel);
        for ii=1:numel(variationaltransversalitycondition.value)
            if ii==1
                if numel(variationaltransversalitycondition.value)>1
                    algebraic.term{ii}=['[' variationaltransversalitycondition.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[variationaltransversalitycondition.value{ii}];
                end
            elseif ii<numel(variationaltransversalitycondition.value)
                algebraic.term{ii}=[variationaltransversalitycondition.value{ii} '; ...'];
            else
                if numel(variationaltransversalitycondition.value)>1
                    algebraic.term{ii}=[variationaltransversalitycondition.value{ii} ']'];
                end
            end
        end
        
    case 'parameterderivativetransversalitycondition'
        parameterderivativetransversalitycondition=retrievemodelinformation(ocStruct,'parameterderivativetransversalitycondition','',getsymkernel);
        for ii=1:numel(parameterderivativetransversalitycondition.value)
            if ii==1
                if numel(parameterderivativetransversalitycondition.value)>1
                    algebraic.term{ii}=['[' parameterderivativetransversalitycondition.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[parameterderivativetransversalitycondition.value{ii}];
                end
            elseif ii<numel(parameterderivativetransversalitycondition.value)
                algebraic.term{ii}=[parameterderivativetransversalitycondition.value{ii} '; ...'];
            else
                if numel(parameterderivativetransversalitycondition.value)>1
                    algebraic.term{ii}=[parameterderivativetransversalitycondition.value{ii} ']'];
                end
            end
        end

    case 'variationalstateconstraint'
        switch classtype
            case 'inequality'
                variationalstateconstraint=retrievemodelinformation(ocStruct,'variationalstateconstraint','',getsymkernel);
                algebraic.term=variationalstateconstraint.value;
            otherwise
        end
        
        
%%%% gradient method
%%%%


    case 'localstatedynamics'
        statedynamics=retrievemodelinformation(ocStruct,'statedynamics','','');
        for ii=1:numel(statedynamics.value)
            if ii==1
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=['[' statedynamics.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[statedynamics.value{ii}];
                end
            elseif ii<numel(statedynamics.value)
                algebraic.term{ii}=[statedynamics.value{ii} '; ...'];
            else
                if numel(statedynamics.value)>1
                    algebraic.term{ii}=[statedynamics.value{ii} ']'];
                end
            end
        end

    case 'localcostatedynamics'
        adjointsystem=retrievemodelinformation(ocStruct,'adjointsystem','','');
        for ii=1:numel(adjointsystem.value)
            if ii==1
                if numel(adjointsystem.value)>1
                    algebraic.term{ii}=['[' adjointsystem.value{ii} '; ...'];
                else
                    algebraic.term{ii}=[adjointsystem.value{ii}];
                end
            elseif ii<numel(adjointsystem.value)
                algebraic.term{ii}=[adjointsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=[adjointsystem.value{ii} ']'];
            end
        end


    case 'pontryaginfunctionDu'
        pontryaginfunctionDu=retrievemodelinformation(ocStruct,'pontryaginfunctionDu');
        numterm=numel(pontryaginfunctionDu.value);
        for ii=1:numterm
            if ii==numterm
                if numterm>1
                    algebraic.term{ii}=['[' pontryaginfunctionDu.value{ii} '; ...'];
                else
                    algebraic.term{ii}=pontryaginfunctionDu.value{ii};
                end
            elseif ii<numterm
                algebraic.term{ii}=[pontryaginfunctionDu.value{ii} '; ...'];
            else
                algebraic.term{ii}=[algebraic.term{ii} '];'];
            end
        end
        
    case 'kuhntuckersolution'
        basename={'lagm','slv','LOCV'};
        kuhntuckersolution=retrievemodelinformation(ocStruct,'kuhntuckersolution','',getsymkernel());
        if ~isempty(kuhntuckersolution.value)
            fn=fieldnames(kuhntuckersolution.value);

            if isempty(fn)
                algebraic.term{1}='[]';
            else
                ctr2=0;
                for ii=1:length(basename)
                    idx=find(strncmp(basename{ii},fn,length(basename{ii})));
                    ctr1=0;
                    for jj=1:length(idx)
                        ctr2=ctr2+1;
                        ctr1=ctr1+1;
                        algebraic.term{ctr2}=[basename{ii} '(:,:,' num2str(ctr1) ')=[']; % kuhntuckersolution.value(1).(fn{idx(jj)}) '; ...'
                        for kk=1:length(kuhntuckersolution.value)
                            ctr2=ctr2+1;
                            if kk<length(kuhntuckersolution.value)
                                algebraic.term{ctr2}=[kuhntuckersolution.value(kk).(fn{idx(jj)}) '; ...'];
                            else
                                algebraic.term{ctr2}=[kuhntuckersolution.value(kk).(fn{idx(jj)}) '];'];
                            end
                        end
                    end
                end
            end
        else
            algebraic.term{1}='out=[];';
        end
        
    case 'specifickuhntuckersolution'
        basename={'lagm','slv','LOCV'};
        kuhntuckersolution=retrievemodelinformation(ocStruct,'specifickuhntuckersolution',arcfield,getsymkernel());
        if ~isempty(kuhntuckersolution.value)
            fn=fieldnames(kuhntuckersolution.value);

            if isempty(fn)
                algebraic.term{1}='[]';
            else
                ctr2=0;
                for ii=1:length(basename)
                    idx=find(strncmp(basename{ii},fn,length(basename{ii})));
                    ctr1=0;
                    for jj=1:length(idx)
                        ctr2=ctr2+1;
                        ctr1=ctr1+1;
                        algebraic.term{ctr2}=[basename{ii} '(:,:,' num2str(ctr1) ')=[']; % kuhntuckersolution.value(1).(fn{idx(jj)}) '; ...'
                        for kk=1:length(kuhntuckersolution.value)
                            ctr2=ctr2+1;
                            if kk<length(kuhntuckersolution.value)
                                algebraic.term{ctr2}=[kuhntuckersolution.value(kk).(fn{idx(jj)}) '; ...'];
                            else
                                algebraic.term{ctr2}=[kuhntuckersolution.value(kk).(fn{idx(jj)}) '];'];
                            end
                        end
                    end
                end
            end
        else
            algebraic.term{1}='out=[];';
        end

    case 'locprojection'
        locprojection=retrievemodelinformation(ocStruct,'locprojection','',getsymkernel());
        if ~isempty(locprojection.value)
            fn=fieldnames(locprojection.value);
            if isempty(fn)
                algebraic.term{1}='[]';
            else
                ctr=0;
                for ii=1:length(fn)
                    for jj=1:length(locprojection.value.(fn{ii}))
                        ctr=ctr+1;
                        algebraic.term{ctr}=locprojection.value.(fn{ii})(jj).explicitexpr;
                    end
                end
            end
        else
            algebraic.term{1}='out=[];';
        end

        % DAE

    case 'canonicalsystemjacobiandae'
        canonicalsystemjacobiandae=retrievemodelinformation(ocStruct,'canonicalsystemjacobiandae','',getsymkernel);
        numterm=numel(canonicalsystemjacobiandae.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemjacobiandae.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemjacobiandae.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemjacobiandae.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemparameterjacobiandae'
        canonicalsystemparameterjacobiandae=retrievemodelinformation(ocStruct,'canonicalsystemparameterjacobiandae','',getsymkernel);
        numterm=numel(canonicalsystemparameterjacobiandae.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[' canonicalsystemparameterjacobiandae.value{ii} '; ...'];
            elseif ii<numterm
                algebraic.term{ii}=[canonicalsystemparameterjacobiandae.value{ii} '; ...'];
            else
                algebraic.term{ii}=canonicalsystemparameterjacobiandae.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'dlagrangefunctiondu'
        dlagrangefunctiondu=retrievemodelinformation(ocStruct,'dlagrangefunctiondu','',getsymkernel);
        numterm=numel(dlagrangefunctiondu.value);
        for ii=1:numterm
            if ii==1
                if ii<numterm
                    algebraic.term{ii}=['[' dlagrangefunctiondu.value{ii} '; ...'];
                else
                    algebraic.term{ii}=dlagrangefunctiondu.value{ii};
                end
            elseif ii<numterm
                algebraic.term{ii}=[dlagrangefunctiondu.value{ii} '; ...'];
            else
                algebraic.term{ii}=dlagrangefunctiondu.value{ii};
            end
        end
        if numterm>1
            algebraic.term{ii}=[algebraic.term{ii} ']'];
        end
    case 'dlagrangefunctiondujacobian'
        dlagrangefunctiondujacobian=retrievemodelinformation(ocStruct,'dlagrangefunctiondujacobian','',getsymkernel);
        numterm=numel(dlagrangefunctiondujacobian.value);
        for ii=1:numterm
            if ii==1
                if ii<numterm
                    algebraic.term{ii}=['[' dlagrangefunctiondujacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['[' dlagrangefunctiondujacobian.value{ii}];
                end
            elseif ii<numterm
                algebraic.term{ii}=[dlagrangefunctiondujacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=dlagrangefunctiondujacobian.value{ii};
            end
        end
            algebraic.term{ii}=[algebraic.term{ii} ']'];


    case 'dlagrangefunctionduparameterjacobian'
        dlagrangefunctionduparameterjacobian=retrievemodelinformation(ocStruct,'dlagrangefunctionduparameterjacobian','',getsymkernel);
        numterm=numel(dlagrangefunctionduparameterjacobian.value);
        for ii=1:numterm
            if ii==1
                if ii<numterm
                    algebraic.term{ii}=['[' dlagrangefunctionduparameterjacobian.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['[' dlagrangefunctionduparameterjacobian.value{ii}];
                end
            elseif ii<numterm
                algebraic.term{ii}=[dlagrangefunctionduparameterjacobian.value{ii} '; ...'];
            else
                algebraic.term{ii}=dlagrangefunctionduparameterjacobian.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'controlconstraintjacobiandae'
        controlconstraintjacobiandae=retrievemodelinformation(ocStruct,'controlconstraintjacobiandae','',getsymkernel);
        numterm=numel(controlconstraintjacobiandae.value);
        for ii=1:numterm
            if ii==1
                if ii<numterm
                    algebraic.term{ii}=['[' controlconstraintjacobiandae.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['[' controlconstraintjacobiandae.value{ii}];
                end
            elseif ii<numterm
                algebraic.term{ii}=[controlconstraintjacobiandae.value{ii} '; ...'];
            else
                algebraic.term{ii}=controlconstraintjacobiandae.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'controlconstraintparameterjacobiandae'
        controlconstraintparameterjacobiandae=retrievemodelinformation(ocStruct,'controlconstraintparameterjacobiandae','',getsymkernel);
        numterm=numel(controlconstraintparameterjacobiandae.value);
        for ii=1:numterm
            if ii==1
                if ii<numterm
                    algebraic.term{ii}=['[' controlconstraintparameterjacobiandae.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['[' controlconstraintparameterjacobiandae.value{ii}];
                end
            elseif ii<numterm
                algebraic.term{ii}=[controlconstraintparameterjacobiandae.value{ii} '; ...'];
            else
                algebraic.term{ii}=controlconstraintparameterjacobiandae.value{ii};
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'objectivefunctionjacobiandae'
        objectivefunctionjacobiandae=retrievemodelinformation(ocStruct,'objectivefunctionjacobiandae','',getsymkernel);
        algebraic.term{1}=['[' objectivefunctionjacobiandae.value{1} '];'];

    case 'objectivefunctionparameterjacobiandae'
        objectivefunctionparameterjacobiandae=retrievemodelinformation(ocStruct,'objectivefunctionparameterjacobiandae','',getsymkernel);
        algebraic.term{1}=['[' objectivefunctionparameterjacobiandae.value{1} '];'];
        
    case 'transversalityconditiondae'
        transversalityconditiondae=retrievemodelinformation(ocStruct,'transversalityconditiondae','',getsymkernel);
        algebraic.term{1}=transversalityconditiondae.value;
end

if vectflag
    if iscell(algebraic.term)
        for ii=1:numel(algebraic.term)
            algebraic.term{ii}=simplevectorize(algebraic.term{ii});
        end
    elseif ischar(algebraic.term)
        algebraic.term=simplevectorize(algebraic.term);
    end
end

function v = simplevectorize(v)
%VECTORIZE Vectorize expression.
%   VECTORIZE(S), when S is a string expression, inserts a '.' before
%   any '^', '*' or '/' in S.  The result is a character string.
%
%   VECTORIZE(FUN), when FUN is an INLINE function object, vectorizes the
%   formula for FUN.  The result is the vectorized version of the INLINE
%   function.
%
%   See also INLINE/FORMULA, INLINE.

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.11 $  $Date: 2002/04/15 04:22:08 $

if isempty(v)
    v = [];
else
    for k = fliplr(find((v=='^') | (v=='*') | (v=='/')))
        v = [v(1:k-1) '.' v(k:end)];
    end
end