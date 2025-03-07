function algebraic=algebraicppdeterm2string(ocStruct,class,vectflag,varargin)
%
% ALGEBRAICPPDETERM2STRING process the algebraic terms of the structure
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

switch class
    case 'spatialcontrolvalue'
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
        spatialcontrolpendenceindex=retrieveppdemodelinformation(ocStruct,'spatialcontroldependenceindex');
        ctr=0;
        for ii=spatialcontrolpendenceindex.value
            ctr=ctr+1;
            if ii==ctr
                if length(spatialcontrolpendenceindex.value)>1
                    algebraic.term{ctr}=['[' optimalvalue.value.(controlname.value{ii}) '; ...'];
                else
                    algebraic.term{ctr}=optimalvalue.value.(controlname.value{ii});
                end
            elseif ii<length(spatialcontrolpendenceindex.value)
                algebraic.term{ctr}=[optimalvalue.value.(controlname.value{ii}) '; ...'];
            else
                algebraic.term{ctr}=[optimalvalue.value.(controlname.value{ii})];
            end
        end
        if length(spatialcontrolpendenceindex.value)>1
            algebraic.term{ctr}=[algebraic.term{ctr} ']'];
        end
        
    case 'nonspatialcontrolvalue'
        gridnumm1=getbasicname('femdatagridnumm1');
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
        spacearg=retrieveppdemodelinformation(ocStruct,'space');
        defaultspacearg=getbasicname('space');
        defaultspaceargmid=getbasicname('spacemid');
        nonspatialcontrolpendenceindex=retrieveppdemodelinformation(ocStruct,'nonspatialcontroldependenceindex');
        ctr=0;
        for ii=nonspatialcontrolpendenceindex.value
            ctr=ctr+1;
            if ii==ctr
                if length(nonspatialcontrolpendenceindex.value)>1
                    nonspatialoptimalvalue=strrep(optimalvalue.value.(controlname.value{ii}),['int_' spacearg.value],'');
                    nonspatialoptimalvalue=strrep(nonspatialoptimalvalue,defaultspacearg,defaultspaceargmid);
                    nonspatialoptimalvalue=['sum(' nonspatialoptimalvalue '*d' defaultspacearg ')'];
                    algebraic.term{ctr}=['[repmat(' nonspatialoptimalvalue ',' gridnumm1 '+1,1); ...'];
                else
                    nonspatialoptimalvalue=strrep(optimalvalue.value.(controlname.value{ii}),['int_' spacearg.value],'');
                    nonspatialoptimalvalue=strrep(nonspatialoptimalvalue,defaultspacearg,defaultspaceargmid);
                    nonspatialoptimalvalue=['sum(' nonspatialoptimalvalue '*d' defaultspacearg ')'];
                    algebraic.term{ctr}=['repmat(' nonspatialoptimalvalue ',' gridnumm1 '+1,1)'];
                end
            elseif ii<length(nonspatialcontrolpendenceindex.value)
                nonspatialoptimalvalue=strrep(optimalvalue.value.(controlname.value{ii}),['int_' spacearg.value],'');
                nonspatialoptimalvalue=strrep(nonspatialoptimalvalue,defaultspacearg,defaultspaceargmid);
                nonspatialoptimalvalue=['sum(' nonspatialoptimalvalue '*d' defaultspacearg ')'];
                algebraic.term{ctr}=['repmat(' nonspatialoptimalvalue ',' gridnumm1 '+1,1); ...'];
            else
                nonspatialoptimalvalue=strrep(optimalvalue.value.(controlname.value{ii}),['int_' spacearg.value],'');
                nonspatialoptimalvalue=strrep(nonspatialoptimalvalue,defaultspacearg,defaultspaceargmid);
                nonspatialoptimalvalue=['sum(' nonspatialoptimalvalue '*d' defaultspacearg ')'];
            end
        end
        if length(nonspatialcontrolpendenceindex.value)>1
            algebraic.term{ctr}=[algebraic.term{ctr} ']'];
        end
        
        
    case 'symboliccontrolvalue'
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        controlname=retrieveppdemodelinformation(ocStruct,'controlname');
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
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalcostatevalue',field2arcidentifier(arcfield),getsymkernel);
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        costatename=retrieveppdemodelinformation(ocStruct,'costatename');
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
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
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
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitycontrolconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
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
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitystateconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrieveppdemodelinformation(ocStruct,'lagrangemultiplierstatename');
        for ii=1:inequalitystateconstraintnum.value
            if ii==1
                if inequalitystateconstraintnum.value>1
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultiplierstatename.value{ii}) '; ...'];
                else
                    algebraic.term{ii}=['[' optimalvalue.value.(lagrangemultiplierstatename.value{ii})];
                end
            elseif ii<inequalitystateconstraintnum.value
                algebraic.term{ii}=[optimalvalue.value.(lagrangemultiplierstatename.value{ii}) '; ...'];
            else
                algebraic.term{ii}=optimalvalue.value.(lagrangemultiplierstatename.value{ii});
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'symboliclagrangemultsc'
        optimalvalue=retrieveppdemodelinformation(ocStruct,'optimalvalue',field2arcidentifier(arcfield));
        inequalitystateconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitystateconstraintnum');
        lagrangemultiplierstatename=retrieveppdemodelinformation(ocStruct,'lagrangemultiplierstatename');
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
        objectiveintegrand=retrieveppdemodelinformation(ocStruct,'objectiveintegrand');
        switch classtype
            case 'current'
                discountfactor=retrieveppdemodelinformation(ocStruct,'discountfactor');
                algebraic.term=[discountfactor.value '*(' objectiveintegrand.value ')'];
            case 'present'
                algebraic.term=objectiveintegrand.value;
        end

    case 'salvagevalue'
        salvagevalue=retrieveppdemodelinformation(ocStruct,'salvagevalue');
        algebraic.term=salvagevalue.value;
        %endtimediscountfactor=retrieveppdemodelinformation(ocStruct,'endtimediscountfactor');
        %algebraic.term=[endtimediscountfactor.value '*(' salvagevalue.value ')'];


    case 'discountedsalvagevalue'
        salvagevalue=retrieveppdemodelinformation(ocStruct,'salvagevalue');
        endtimediscountfactor=retrieveppdemodelinformation(ocStruct,'endtimediscountfactor');
        if strcmp(salvagevalue.value,'0')
            algebraic.term=salvagevalue.value;
        else
            algebraic.term=[endtimediscountfactor.value '*(' salvagevalue.value ')'];
        end
    case 'controlconstraint'
        switch classtype
            case 'inequality'
                inequalitycontrolconstraint=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraint');
                algebraic.term=inequalitycontrolconstraint.value;
            otherwise
        end

    case 'constraint'
        inequalitycontrolconstraint=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraint');
        inequalitystateconstraint=retrieveppdemodelinformation(ocStruct,'inequalitystateconstraint');
        algebraic.term=[inequalitycontrolconstraint.value inequalitystateconstraint.value];
        
    case 'canonicalsystem'
        isexplicitspatial=retrieveppdemodelinformation(ocStruct,'isexplicitspatial');
        if isexplicitspatial.value
            space=retrieveppdemodelinformation(ocStruct,'space');
            dependent=getbasicname('dependent');
        end
        canonicalsystem=retrieveppdemodelinformation(ocStruct,'canonicalsystem','','');
        for ii=1:numel(canonicalsystem.value)
            if isexplicitspatial.value
                cansysval=regexprep(canonicalsystem.value{ii},['\<' space.value '\>'],['repmat(' space.value ',1,size(' dependent ',2))']);
            else
                cansysval=canonicalsystem.value{ii};
            end
            if ii==1
                algebraic.term{ii}=['[' cansysval '; ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=[cansysval '; ...'];
            else
                algebraic.term{ii}=[cansysval ']'];
            end
        end

    case 'canonicalsystemcont'
        canonicalsystem=retrieveppdemodelinformation(ocStruct,'canonicalsystemcont','','');
        for ii=1:numel(canonicalsystem.value)
            if ii==1
                algebraic.term{ii}=['[' canonicalsystem.value{ii} '; ...'];
            elseif ii<numel(canonicalsystem.value)
                algebraic.term{ii}=[canonicalsystem.value{ii} '; ...'];
            else
                algebraic.term{ii}=[canonicalsystem.value{ii} ']'];
            end
        end

    case 'exogenousfunctionterm'
        exogenousfunctionterm=retrieveppdemodelinformation(ocStruct,'exogenousfunctionterm');
        exogenousfunctionnum=retrieveppdemodelinformation(ocStruct,'exogenousfunctionnum');
        for ii=1:exogenousfunctionnum.value
            if ii==1
                if exogenousfunctionnum.value>1
                    algebraic.term{ii}=['out=[' exogenousfunctionterm.value{ii} '; ...'];
                else
                    algebraic.term{ii}=['out=' exogenousfunctionterm.value{ii} ';'];
                end
            elseif ii<exogenousfunctionnum.value
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} '; ...'];
            else
                algebraic.term{ii}=[exogenousfunctionterm.value{ii} ']'];
            end
        end

    case 'algebraicequation'
        arcidentifier=field2arcidentifier(arcfield);
        algebraicequation=retrieveppdemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrieveppdemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
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
        algebraicequation=retrieveppdemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrieveppdemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
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
        algebraicequationimplicit=retrieveppdemodelinformation(ocStruct,'algebraicequationimplicit',arcidentifier,getsymkernel);
        algebraicequationimplicitnum=numel(algebraicequationimplicit.value);
        for ii=1:algebraicequationimplicitnum
            if ii==1
                if algebraicequationimplicitnum>1
                    algebraic.term{ii}=['[' algebraicequationimplicit.value{ii} '; ...'];
                else
                    algebraic.term{ii}=algebraicequationimplicit.value{ii};
                end
            elseif ii<numel(algebraicequationimplicit.value)
                algebraic.term{ii}=['' algebraicequationimplicit.value{ii} '; ...'];
            else
                algebraic.term{ii}=['' algebraicequationimplicit.value{ii} ']'];
            end
        end

    case 'symboliccanonicalsystem'
        canonicalsystem=retrieveppdemodelinformation(ocStruct,'canonicalsystem','','');
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
        algebraicequation=retrieveppdemodelinformation(ocStruct,'algebraicequation',arcidentifier,getsymkernel);
        algebraicequationnum=retrieveppdemodelinformation(ocStruct,'algebraicequationnum',arcidentifier,'');
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
        canonicalsystemjacobian=retrieveppdemodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
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
        canonicalsystemjacobian=retrieveppdemodelinformation(ocStruct,'canonicalsystemjacobianstatecontrol',arcidentifier,getsymkernel);
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
        specificcanonicalsysteminstatecontrol=retrieveppdemodelinformation(ocStruct,'specificcanonicalsysteminstatecontrol',arcidentifier,getsymkernel);
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
        specificcanonicalsystem=retrieveppdemodelinformation(ocStruct,'specificcanonicalsystem',arcidentifier,getsymkernel);
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
        specificpontryaginfunction=retrieveppdemodelinformation(ocStruct,'specificpontryaginfunction',arcidentifier,getsymkernel);
        
        algebraic.term=specificpontryaginfunction.value;

    case 'symboliccanonicalsystem'
        arcidentifier=field2arcidentifier(arcfield);
        specificcanonicalsystem=retrieveppdemodelinformation(ocStruct,'symboliccanonicalsystem',arcidentifier,getsymkernel);
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
        canonicalsystemjacobian=retrieveppdemodelinformation(ocStruct,'statecostatejacobian',arcidentifier,getsymkernel);
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
        gridnum=getbasicname('femdatagridnum');
        canonicalsystemjacobian=retrieveppdemodelinformation(ocStruct,'canonicalsystemjacobian',arcidentifier,getsymkernel);
        numterm=numel(canonicalsystemjacobian.value);
        for ii=1:numterm
            if ii==1
                canonicalsystemjacobianvalue=regexp(canonicalsystemjacobian.value{ii},',','split');
                canonicalsystemjacobianline=[];
                for jj=1:length(canonicalsystemjacobianvalue)
                    canonicalsystemjacobianvalue{jj}=int2sum(canonicalsystemjacobianvalue{jj},ocStruct);

                    canonicalsystemjacobianline=[canonicalsystemjacobianline 'spdiags(' canonicalsystemjacobianvalue{jj} ',0,' gridnum ',' gridnum '),'];
                end
                canonicalsystemjacobianline(end)=[];
                algebraic.term{ii}=['[' canonicalsystemjacobianline '; ...'];
            elseif ii<numterm
                canonicalsystemjacobianvalue=regexp(canonicalsystemjacobian.value{ii},',','split');
                canonicalsystemjacobianline=[];
                for jj=1:length(canonicalsystemjacobianvalue)
                    canonicalsystemjacobianvalue{jj}=int2sum(canonicalsystemjacobianvalue{jj},ocStruct);                    canonicalsystemjacobianline=[canonicalsystemjacobianline 'spdiags(' canonicalsystemjacobianvalue{jj} ',0,' gridnum ',' gridnum '),'];
                end
                canonicalsystemjacobianline(end)=[];
                algebraic.term{ii}=[canonicalsystemjacobianline '; ...'];
            else
                canonicalsystemjacobianvalue=regexp(canonicalsystemjacobian.value{ii},',','split');
                canonicalsystemjacobianline=[];
                for jj=1:length(canonicalsystemjacobianvalue)
                    canonicalsystemjacobianvalue{jj}=int2sum(canonicalsystemjacobianvalue{jj},ocStruct);                    canonicalsystemjacobianline=[canonicalsystemjacobianline 'spdiags(' canonicalsystemjacobianvalue{jj} ',0,' gridnum ',' gridnum '),'];
                end
                canonicalsystemjacobianline(end)=[];
                algebraic.term{ii}=canonicalsystemjacobianline;
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'canonicalsystemjacobianstatecontrol'
        arcidentifier=field2arcidentifier(arcfield);
        canonicalsystemjacobianstatecontrol=retrieveppdemodelinformation(ocStruct,'canonicalsystemjacobianstatecontrol',arcidentifier,getsymkernel);
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
        optimalcontroldynamicsleftside=retrieveppdemodelinformation(ocStruct,'optimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'totalimplicitvariableindex');
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
        optimalcontroldynamicsleftside=retrieveppdemodelinformation(ocStruct,'optimalcontroldynamics',arcidentifier,getsymkernel);
        %implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'totalimplicitvariableindex',arcidentifier);
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
        optimalcontroldynamicsleftside=retrieveppdemodelinformation(ocStruct,'optimalcontroldynamics',arcidentifier,getsymkernel);
        %implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
        implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'totalimplicitvariableindex',arcidentifier);
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
        optimalcontroldynamicsleftside=retrieveppdemodelinformation(ocStruct,'invoptimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
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
        tensoroptimalcontroldynamicsleftside=retrieveppdemodelinformation(ocStruct,'tensoroptimalcontroldynamicsleftside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
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
        jacobianoptimalcontroldynamicsrightside=retrieveppdemodelinformation(ocStruct,'jacobianoptimalcontroldynamicsrightside',arcidentifier,getsymkernel);
        implicitnonlinearcontrolindex=retrieveppdemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier);
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
        pontryaginfunctionDuDX=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDuDX','','');
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
        pontryaginfunctionDlmmcDX=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDlmmcDX',arcidentifier,getsymkernel);
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
        canonicalsystemparameterjacobian=retrieveppdemodelinformation(ocStruct,'canonicalsystemparameterjacobian',arcidentifier,getsymkernel);
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
        statecontrolcanonicalsystemparameterjacobian=retrieveppdemodelinformation(ocStruct,'statecontrolcanonicalsystemparameterjacobian',arcidentifier,getsymkernel);
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
        canonicalsystemderivativetime=retrieveppdemodelinformation(ocStruct,'canonicalsystemderivativetime',arcidentifier,getsymkernel);
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
        canonicalsystemhessian=retrieveppdemodelinformation(ocStruct,'canonicalsystemhessian',arcidentifier,getsymkernel);
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
        canonicalsystemparameterhessian=retrieveppdemodelinformation(ocStruct,'canonicalsystemtotalhessian',arcidentifier,getsymkernel);
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

    case 'parametervariables'
        parametername=retrieveppdemodelinformation(ocStruct,'parametername');
        algebraic.term=parametername.value;

    case 'discobjectivefunction'
        discountfactor=retrieveppdemodelinformation(ocStruct,'discountfactor');
        gridnum=getbasicname('femdatagridnum');
        
        discobjectivefunction=retrieveppdemodelinformation(ocStruct,'discobjectivefunction');
        discobjectivefunction.value=strrep(discobjectivefunction.value,discountfactor.value,['repmat(' discountfactor.value ',' gridnum ',1)']);
        algebraic.term{1}=discobjectivefunction.value;

    case 'discobjectivefunctionderivativetime'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionderivativetime=retrieveppdemodelinformation(ocStruct,'discobjectivefunctionderivativetime',arcidentifier,getsymkernel);
        algebraic.term=discobjectivefunctionderivativetime.value;


    case 'symbolicdiscobjectivefunction'
        discobjectivefunction=retrieveppdemodelinformation(ocStruct,'discobjectivefunction');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        algebraic.term{1}=[symstr '(''' discobjectivefunction.value ''')'];

    case 'discobjectivefunctionjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrieveppdemodelinformation(ocStruct,'discobjectivefunctionDX',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        gridnum=getbasicname('femdatagridnum');
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[spdiags(' discobjectivefunctionjacobian.value{ii} ',0,' gridnum ',' gridnum '), ...'];
            elseif ii<numterm
                algebraic.term{ii}=['spdiags(' discobjectivefunctionjacobian.value{ii} ',0,' gridnum ',' gridnum '), ...'];
            else
                algebraic.term{ii}=['spdiags(' discobjectivefunctionjacobian.value{ii} ',0,' gridnum ',' gridnum ')'];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'discobjectivefunctionparameterjacobian'
        arcidentifier=field2arcidentifier(arcfield);
        discobjectivefunctionjacobian=retrieveppdemodelinformation(ocStruct,'discobjectivefunctionDP',arcidentifier,getsymkernel);
        numterm=numel(discobjectivefunctionjacobian.value);
        for ii=1:numterm
            if ii==1
                algebraic.term{ii}=['[spdiags(' discobjectivefunctionjacobian.value{ii} '), ...'];
            elseif ii<numterm
                algebraic.term{ii}=['spdiags(' discobjectivefunctionjacobian.value{ii} '), ...'];
            else
                algebraic.term{ii}=['spdiags(' discobjectivefunctionjacobian.value{ii} ')'];
            end
        end
        algebraic.term{ii}=[algebraic.term{ii} ']'];

    case 'hamiltonianfunction'
        hamiltonianfunction=retrieveppdemodelinformation(ocStruct,'hamiltonianfunction');
        algebraic.term{1}=hamiltonianfunction.value;

    case 'pontryaginfunction'
        pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
        algebraic.term{1}=pontryaginfunction.value;

    case 'symbolicpontryaginfunction'
        pontryaginfunction=retrieveppdemodelinformation(ocStruct,'pontryaginfunction');
        if verLessThan('symbolic','8')
            symstr='mystr2sym';
        else
            symstr='mystr2sym';
        end
        algebraic.term{1}=[symstr '(''' pontryaginfunction.value ''')'];

    case 'pontryaginfunctionDx'
        pontryaginfunctionDx=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDx');
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

    case 'pontryaginfunctionDu2'
        pontryaginfunctionDu2=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDu2');
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

    case 'pontryaginfunctionDX'
        arcidentifier=field2arcidentifier(arcfield);
        pontryaginfunctionDx=retrieveppdemodelinformation(ocStruct,'pontryaginfunctionDX',arcidentifier,getsymkernel);
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

    case 'inequalitystateconstrainttimederivative'
        inequalitystateconstraint=retrieveppdemodelinformation(ocStruct,'inequalitystateconstraint');
        inequalitystateconstrainttimederivative=retrieveppdemodelinformation(ocStruct,'inequalitystateconstrainttimederivative');
        inequalitystateconstraintorder=retrieveppdemodelinformation(ocStruct,'inequalitystateconstraintorder');
        inequalitystateconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitystateconstraintnum');
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
        hamiltonianDu=retrieveppdemodelinformation(ocStruct,'hamiltonianDu','',getsymkernel);
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
        statedynamics=retrieveppdemodelinformation(ocStruct,'statedynamics','','');
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
        adjointsystem=retrieveppdemodelinformation(ocStruct,'adjointsystem','0','');
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
        objectivefunction=retrieveppdemodelinformation(ocStruct,'objectiveintegrand');
        isexplicitspatial=retrieveppdemodelinformation(ocStruct,'isexplicitspatial');
        if isexplicitspatial.value
            space=retrieveppdemodelinformation(ocStruct,'space');
            gridnum=getbasicname('femdatagridnum');
        end
        if isexplicitspatial.value
            algebraic.term{1}=regexprep(objectivefunction.value,['\<' space.value '\>'],['repmat(' space.value ',' gridnum ',1)']);
        else
            algebraic.term{1}=objectivefunction.value;
        end
        
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

function str=int2sum(str,ocStruct)

gridnumm1=getbasicname('femdatagridnumm1');
spacearg=retrieveppdemodelinformation(ocStruct,'space');
defaultspacearg=getbasicname('space');
defaultmidspacearg=getbasicname('spacemid');

removearg={['D(int_' spacearg.value ')'],['int_' spacearg.value]};

for ii=1:length(removearg)
    idx4removearg=findstr(str,removearg{ii});
    if ~isempty(idx4removearg)
        idx4removearg=unique([idx4removearg length(str)]);
        lengthremovearg=length(removearg{ii});
        out='';
        initidx=1;
        for jj=1:length(idx4removearg)-1
            actualstring=str(idx4removearg(jj)+lengthremovearg:idx4removearg(jj+1));
            paranthesisloc=zeros(1,length(actualstring));
            paranthesisloc(findstr(actualstring,'('))=1;
            paranthesisloc(findstr(actualstring,')'))=-1;
            paranthesisloc=cumsum(paranthesisloc);
            closingparidx=find(~paranthesisloc,1);
            tmpstring=strrep(actualstring(1:closingparidx),defaultspacearg,defaultmidspacearg);
            tmpstring=['repmat(sum(' tmpstring '*d' defaultspacearg '),' gridnumm1 '+1,1)'];
            out=[out str(initidx:idx4removearg(jj)-1) tmpstring];
            initidx=idx4removearg(jj)+lengthremovearg+closingparidx;
        end
        out=[out str(initidx:end)];
        str=out;
    end
end