function varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,varname,overwrite,varargin)
%
% a structure is generated which is then used to replace variables in
% template model files
% with overwrite 1 an already existing field of varStruct is overwritten
% the structure consists of a field given by the name of the variable. If
% the actual value of the variable depends on the arc the field string
% contains a cell array of strings and the field arcependent is set to one.
% multline 1 ... lines>1 are intended, e.g., CANONICALSYSTEMDYNAMICS
% multline 2 ...lines are not intended, e.g., PARAMETERVALUES
opt=[];
if nargin<=3
    overwrite=0;
end
if ~overwrite
    if isfield(varStruct,varname)
        return
    end
end
if nargin>=5
    opt=varargin{1};
end

varStruct.(varname).arcdependent=0;
varStruct.(varname).arcidentifier=[];
varStruct.(varname).vectorize=0;
varStruct.(varname).type='';
varStruct.(varname).multline=0;

switch varname
    % option dependent variables
    case 'EXPLICITCANONICALSYSTEM'
        if ~isempty(opt) && strcmp(getocoptions(opt,'INIT','ControlDynamics'),'implicit')
            varStruct.(varname).string='0';
        else
            varStruct.(varname).string='1';
        end

        % model dependent variables
    case 'INDEPENDENT'
        varStruct.(varname).string=getbasicname('independent');

    case 'DEPENDENTVAR'
        varStruct.(varname).string=getbasicname('dependent');

    case 'ENDTIME'
        varStruct.(varname).string=getbasicname('endtime');

    case 'PARVAR'
        varStruct.(varname).string=getbasicname('parametervariables');

    case 'ARCVAR'
        varStruct.(varname).string=getbasicname('arcidentifiervar');

    case 'SUBSFLAG'
        varStruct.(varname).string='s';

    case 'LAGRANGEMULTCC'
        varStruct.(varname).string=getbasicname('lagrangemultcc');

    case 'LAGRANGEMULTSC'
        varStruct.(varname).string=getbasicname('lagrangemultsc');

    case 'CONTROL'
        varStruct.(varname).string=getbasicname('control');

    case 'CONSTRAINTVARIABLE'
        varStruct.(varname).string=getbasicname('constraint');
        
    case 'CONTROLDYNAMICSVAR'
        varStruct.(varname).string=['D' getbasicname('control') 'Dt'];

    case 'STATE'
        varStruct.(varname).string=getbasicname('state');

    case 'COSTATE'
        varStruct.(varname).string=getbasicname('costate');

    case 'JUMPARG'
        varStruct.(varname).string='jumparg';

    case 'VJUMPARG'
        %variational jumparg
        varStruct.(varname).string='vjumparg';

    case 'JUMPID'
        varStruct.(varname).string='jumpid';
        
    case 'EXOGENOUSFUNCTION'
        varStruct.(varname).string=getbasicname('exogenousfunction');

    case 'EXOGENOUSFUNCTIONDX'
        varStruct.(varname).string=['D1' getbasicname('exogenousfunction') '_Dx'];

    case 'EXOGENOUSFUNCTIONNUM'
        exogenousfunctionnum=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
        varStruct.(varname).string=num2str(exogenousfunctionnum.value);
        

    case 'EXOGENOUSDYNAMICSTERM'
        try
            exogenousfunctionterm=algebraicterm2string(ocStruct,'exogenousdynamicsterm',1,'');
            for jj=1:numel(exogenousfunctionterm.term)
                if jj==1
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                else
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSDYNAMICSJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousdynamicsjacobianterm=algebraicterm2string(ocStruct,'exogenousdynamicsjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousdynamicsjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousdynamicsjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousdynamicsjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSSDYNAMICSPARAMETERJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousdynamicsparameterjacobianterm=algebraicterm2string(ocStruct,'exogenousdynamicsparameterjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousdynamicsparameterjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousdynamicsparameterjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousdynamicsparameterjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSFUNCTIONTERM'
        try
            exogenousfunctionterm=algebraicterm2string(ocStruct,'exogenousfunctionterm',1,'');
            for jj=1:numel(exogenousfunctionterm.term)
                if jj==1
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                else
                    varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousjacobianterm=algebraicterm2string(ocStruct,'exogenousjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSPARAMETERJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousparameterjacobianterm=algebraicterm2string(ocStruct,'exogenousparameterjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousparameterjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousparameterjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousparameterjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSCTJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousctjacobianterm=algebraicterm2string(ocStruct,'exogenousctjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousctjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousctjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousctjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSDYNAMICS'
        varStruct.(varname).string=getbasicname('exogenousdynamics');

    case 'EXOGENOUSDYNAMICSDX'
        varStruct.(varname).string=['D1' getbasicname('exogenousdynamics') '_Dx'];

    case 'EXOGENOUSDYNAMICSNUM'
        exogenousdynamicsnum=retrievemodelinformation(ocStruct,'exogenousdynamicsnum');
        varStruct.(varname).string=num2str(exogenousdynamicsnum.value);

    case 'EXOGENOUSDYNAMICSTERM'
        try
            exogenousdynamicsterm=algebraicterm2string(ocStruct,'exogenousdynamicsterm',1,'');
            for jj=1:numel(exogenousdynamicsterm.term)
                if jj==1
                    varStruct.(varname).string{jj}=exogenousdynamicsterm.term{jj};
                else
                    varStruct.(varname).string{jj}=exogenousdynamicsterm.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousjacobianterm=algebraicterm2string(ocStruct,'exogenousjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        

    case 'EXOGENOUSPARAMETERJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousparameterjacobianterm=algebraicterm2string(ocStruct,'exogenousparameterjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousparameterjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousparameterjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousparameterjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXOGENOUSCTJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                exogenousctjacobianterm=algebraicterm2string(ocStruct,'exogenousctjacobianterm',0,arcfield{ii});
                for jj=1:numel(exogenousctjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' exogenousctjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=exogenousctjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ADMISSIBLEFACTOR'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                admissiblefactor=algebraicterm2string(ocStruct,'admissiblefactor',0,arcfield{ii});
                varStruct.(varname).string{ii}=['out=' admissiblefactor.term ';'];
            end
            varStruct.(varname).multline=0;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'NONSMOOTHFUNCTION'
        varStruct.(varname).string=getbasicname('nonsmoothfunction');
        
    case 'NONSMOOTHFUNCTIONDX'
        varStruct.(varname).string=['D1' getbasicname('nonsmoothfunction') '_Dx'];

    case 'NONSMOOTHFUNCTIONNUM'
        nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
        varStruct.(varname).string=num2str(nonsmoothfunctionnum.value);

    case 'SWITCHINGVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                switchingvalue=algebraicterm2string(ocStruct,'switchingvalue',1,arcfield{ii});
                for jj=1:numel(switchingvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' switchingvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=switchingvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'GENERATIONDATE'
        varStruct.(varname).string=datestr(now);

    case 'MODELNAME'
        varStruct.(varname).string=modelname(ocStruct);

    case 'UPPERMODELNAME'
        varStruct.(varname).string=upper(modelname(ocStruct));

    case 'ARCIDENTIFIER'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        varStruct.(varname).string=arcidentifier.value;

    case 'NUMBEROFARCS'
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        varStruct.(varname).string=num2str(arcnum.value);

    case 'ARCARGUMENT'
        arcargument=retrievemodelinformation(ocStruct,'argument');
        varStruct.(varname).string=mat2str(arcargument.value);

    case 'ARCCOLOR'
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        colour=get(0,'DefaultAxesColorOrder');
        if arcnum.value<=size(colour,1)
            varStruct.(varname).string=mat2str(colour(1:arcnum.value,:));
        else
            varStruct.(varname).string=mat2str([colour;rand(arcnum.value-size(colour,1),3)]);
        end

    case 'EDGES'
        arcargument=retrievemodelinformation(ocStruct,'argument');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        if arcnum.value>1
            nc=nchoosek(arcargument.value,2);
            varStruct.(varname).string=mat2str([nc;nc(:,[2 1])]);
        else
            varStruct.(varname).string='[]';
        end

    case 'CONTROLVALUESUBSTITUTION'
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        indexstr='';
        for ii=1:controlnum.value
            indexstr=[indexstr '''' controlname.value{ii} ''','];
        end
        coordstr=basename2vectorstring(getbasicname('control'),1:controlnum.value,'coord');
        if inequalitycontrolconstraintnum.value
            for ii=1:inequalitycontrolconstraintnum.value
                indexstr=[indexstr '''' lagrangemultipliercontrolname.value{ii} ''','];
            end
            %indexstr=[indexstr ''',''' basename2vectorstring(getbasicname('lagrangemultcc'),1:inequalitycontrolconstraintnum.value,'index',''',''')];
            coordstr=[coordstr ',' basename2vectorstring(getbasicname('lagrangemultcc'),1:inequalitycontrolconstraintnum.value,'coord')];
        end
        if ~isempty(indexstr)
            indexstr(end)=[];
            varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];
        else
            varStruct.(varname).string=indexstr;
        end

    case 'LAGRANGEMULTIPLIERSUBSTITUTION'
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        lmname=getbasicname('lagrangemultcc');
        if inequalitycontrolconstraintnum.value
            str=sprintf('\tfor ii=1:%d\n\t\tif %s(ii)==0\n\t\t\tout=subs(out,[''%s'' num2str(ii)],0);\n\t\tend\n\tend', ...
                inequalitycontrolconstraintnum.value,lmname,lmname);
        else
            str='';
        end
        varStruct.(varname).string=str;
    case 'COSTATEVALUESUBSTITUTION'
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        indexstr='';
        for ii=1:costatenum.value
            indexstr=[indexstr '''' costatename.value{ii} ''','];
        end
        coordstr=basename2vectorstring(getbasicname('costate'),1:costatenum.value,'coord');
        indexstr(end)=[];
        varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];

    case 'STATEVALUESUBSTITUTION'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        statename=retrievemodelinformation(ocStruct,'statename');
        indexstr='';
        for ii=1:statenum.value
            indexstr=[indexstr '''' statename.value{ii} ''','];
        end
        coordstr=basename2vectorstring(getbasicname('state'),1:statenum.value,'coord');
        indexstr(end)=[];
        varStruct.(varname).string=['out=subs(out,{' indexstr '},{' coordstr '});'];

    case 'CONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value+inequalitystateconstraintnum.value);

    case 'CONTROLCONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value);

    case 'MIXEDCONTROLCONSTRAINTNUM'
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitymixedcontrolconstraintnum');
        varStruct.(varname).string=num2str(inequalitycontrolconstraintnum.value);

    case 'ARCWITHSTATECONSTRAINT'
        arcwithinequalitystateconstraint=retrievemodelinformation(ocStruct,'arcwithinequalitystateconstraint');
        varStruct.(varname).string=arcwithinequalitystateconstraint.value;

    case 'STATECONSTRAINTNUM'
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        varStruct.(varname).string='';
        varStruct.(varname).string=num2str(inequalitystateconstraintnum.value);

    case 'STATECONSTRAINTORDER'
        inequalitystateconstraintorder=retrievemodelinformation(ocStruct,'inequalitystateconstraintorder');
        try
            for ii=0:max(inequalitystateconstraintorder.value)
                varStruct.(varname).string{ii+1}=num2str(ii);
            end
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=[];
        end
       
    case 'STATECONSTRAINTBCR'
        %inequalitystateconstraintorder=retrievemodelinformation(ocStruct,'inequalitystateconstraintorder');
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        time=getbasicname('independent');
        dependent=getbasicname('dependent');
        parametervariables=getbasicname('parametervariables');
        for ii=1:arcnum.value
            fidx=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier.value{ii});
            if ~isempty(fidx.value)
                varStruct.(varname).string{ii}{1}=['DscDtN=' modelname(ocStruct) 'StateConstraintTotalTimeDerivative(' time ',' dependent 'r,' parametervariables ',edge(2),1);'];
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}{2}=['out=DscDtN(' num2str(fidx.value) ');'];
                else
                    varStruct.(varname).string{ii}{2}=['out=DscDtN([' num2str(fidx.value) ']);'];
                end
            else
                varStruct.(varname).string{ii}{1}='out=[];';
            end
        end
        varStruct.(varname).arcdependent=1;
        varStruct.(varname).multline=1;

    case 'STATECONSTRAINTBCL'
        %inequalitystateconstraintorder=retrievemodelinformation(ocStruct,'inequalitystateconstraintorder');
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        time=getbasicname('independent');
        dependent=getbasicname('dependent');
        parametervariables=getbasicname('parametervariables');
        for ii=1:arcnum.value
            fidx=retrievemodelinformation(ocStruct,'nonzerolmscindex',arcidentifier.value{ii});
            if ~isempty(fidx.value)
                varStruct.(varname).string{ii}{1}=['DscDtN=' modelname(ocStruct) 'StateConstraintTotalTimeDerivative(' time ',' dependent 'l,' parametervariables ',edge(1),0);'];
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}{2}=['out=[out; DscDtN(' num2str(fidx.value) ')];'];
                else
                    varStruct.(varname).string{ii}{2}=['out=[out; DscDtN([' num2str(fidx.value) '])];'];
                end
                varStruct.(varname).string{ii}{3}=['DscDtN=' modelname(ocStruct) 'StateConstraintTotalTimeDerivative(' time ',' dependent 'l,' parametervariables ',edge(1),1);'];
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}{4}=['out=[out; DscDtN(' num2str(fidx.value) ')];'];
                else
                    varStruct.(varname).string{ii}{4}=['out=[out; DscDtN([' num2str(fidx.value) '])];'];
                end
            end
        end
        varStruct.(varname).arcdependent=1;
        varStruct.(varname).multline=1;

    case 'STATECONSTRAINTTOTALTIMEDERIVATIVE'
        inequalitystateconstrainttimederivative=algebraicterm2string(ocStruct,'inequalitystateconstrainttimederivative',1,'');
        inequalitystateconstraintnum=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
        try
            for ii=1:numel(inequalitystateconstrainttimederivative.term)
                if inequalitystateconstraintnum.value==1
                    varStruct.(varname).string{ii}=['out=' inequalitystateconstrainttimederivative.term{ii}];
                else
                    varStruct.(varname).string{ii}=['out=[' inequalitystateconstrainttimederivative.term{ii} ']'];
                end
                varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            end
            %end
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ISAUTONOMOUS'
        autonomousflag=retrievemodelinformation(ocStruct,'autonomous');
        varStruct.(varname).string=num2str(autonomousflag.value);

    case 'IMPLICITCONTROLS'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);

 
    case 'NUMERICJACOBIAN'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);
 
    case 'NUMERICHESSIAN'
        implicitcontrolsflag=implicitcontrols(ocStruct);
        varStruct.(varname).string=num2str(implicitcontrolsflag);
        
    case 'CONSTRAINT'
        try
            constraint=algebraicterm2string(ocStruct,'constraint',1);
            varStruct.(varname).string{1}='out=[';
            for jj=1:numel(constraint.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=[' constraint.term{jj}];
                    if numel(constraint.term)>1
                        varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '; ...'];
                    end
                elseif jj<numel(constraint.term)
                    varStruct.(varname).string{jj}=[constraint.term{jj} '; ...'];
                else
                    varStruct.(varname).string{jj}=constraint.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '];'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCONSTRAINTTERM'
        try
            sybmolicconstraint=algebraicterm2string(ocStruct,'symbolicconstraint',0,'');
            for jj=1:numel(sybmolicconstraint.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' sybmolicconstraint.term{jj}];
                else
                    varStruct.(varname).string{jj}=sybmolicconstraint.term{jj};
                end
            end
            if ~isempty(jj)
                varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            else
                varStruct.(varname).string{1}='out=[];';
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECONSTRAINT'
        try
            inequalitystateconstraint=algebraicterm2string(ocStruct,'stateconstraint',1,[],'inequality');
            for jj=1:numel(inequalitystateconstraint.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=[' inequalitystateconstraint.term{jj}];
                    if numel(inequalitystateconstraint.term)>1
                        varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '; ...'];
                    end
                elseif jj<numel(inequalitystateconstraint.term)
                    varStruct.(varname).string{jj}=[inequalitystateconstraint.term{jj} '; ...'];
                else
                    varStruct.(varname).string{jj}=inequalitystateconstraint.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '];'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECOSTATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=['1:' int2str(statenum.value+costatenum.value)];

    case 'CANONICALSYSTEMEQUATIONCOORD'
        canonicalsystemequationnum=retrievemodelinformation(ocStruct,'canonicalsystemequationnum');
        varStruct.(varname).string=['1:' int2str(canonicalsystemequationnum.value)];

    case 'IMPLICITCONTROLNUM'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            equationnum=zeros(1,numel(arcfield));
            counter=1;
            while counter<=numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{counter});
                implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier{1});
                equationnum(counter)=length(implicitnonlinearcontrolindex.value);
                counter=counter+1;
            end
            varStruct.(varname).string=sprintf('%d,',equationnum);
            varStruct.(varname).string(end)=[];
            varStruct.(varname).string=['[' varStruct.(varname).string ']'];
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'GENERALIZEDCONTROLNUM'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            equationnum=zeros(1,numel(arcfield));
            counter=1;
            while counter<=numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{counter});
                generalizedcontrolnum=retrievemodelinformation(ocStruct,'generalizedcontrolnum',arcidentifier{1});
                equationnum(counter)=generalizedcontrolnum.value;
                counter=counter+1;
            end
            varStruct.(varname).string=sprintf('%d,',equationnum);
            varStruct.(varname).string(end)=[];
            varStruct.(varname).string=['[' varStruct.(varname).string ']'];
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONTROLNUM'
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        varStruct.(varname).string=int2str(controlnum.value);

    case 'STATECOSTATENUM'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=int2str(statenum.value+costatenum.value);

    case 'STATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        varStruct.(varname).string=['1:' int2str(statenum.value)];

    case 'STATENUM'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        varStruct.(varname).string=int2str(statenum.value);

    case 'PARAMETERNUM'
        parameternum=retrievemodelinformation(ocStruct,'parameternum');
        varStruct.(varname).string=int2str(parameternum.value);

    case 'DEFAULTPLOTCOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string=['h=plot(' depvarname '(1,:),' depvarname '(2,:));'];

    case 'DEFAULTPLOTTHREECOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string=['h=plot3(' depvarname '(1,:),' depvarname '(2,:),' depvarname '(3,:));'];

    case 'DEFAULTPLOTINDIFFERENCECOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string{1}=['h=plot(' depvarname '(1,leftarcindex(1):rightarcindex(1)),' depvarname '(2,leftarcindex(1):rightarcindex(1)), ...'];
        varStruct.(varname).string{2}=['    ' depvarname '(1,leftarcindex(2):rightarcindex(2)),' depvarname '(2,leftarcindex(2):rightarcindex(2)));'];
        varStruct.(varname).multline=1;

    case 'DEFAULTPLOTTHREEINDIFFERENCECOMMAND'
        depvarname=getbasicname('dependent');
        varStruct.(varname).string{1}=['h=plot3(' depvarname '(1,leftarcindex(1):rightarcindex(1)),' depvarname '(2,leftarcindex(1):rightarcindex(1)),' depvarname '(3,leftarcindex(1):rightarcindex(1)), ...'];
        varStruct.(varname).string{2}=['    ' depvarname '(1,leftarcindex(2):rightarcindex(2)),' depvarname '(2,leftarcindex(2):rightarcindex(2)),' depvarname '(3,leftarcindex(2):rightarcindex(2)));'];
        varStruct.(varname).multline=1;

    case 'LATEXUSERDEPENDENTNAMEFIRST'
        statename=retrievemodelinformation(ocStruct,'statename');
        [val1 val2 val3]=regexp(statename.value{1},'([0-9]+)\>','split');
        if verLessThan('symbolic','8')
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        else
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        end
        %varStruct.(varname).string=['$' regexprep(statename.value{1},'([0-9]+)\>','_$1') '$'];

    case 'LATEXUSERDEPENDENTNAMESECOND'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        if statenum.value>1
            statename=retrievemodelinformation(ocStruct,'statename');
            labelname=statename.value{2};
        else
            costatename=retrievemodelinformation(ocStruct,'costatename');
            labelname=costatename.value{1};
        end
        [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
        if verLessThan('symbolic','8')
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        else
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        end

    case 'LATEXUSERDEPENDENTNAMETHIRD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        if statenum.value<=2
            labelname='';
        else
            statename=retrievemodelinformation(ocStruct,'statename');
            labelname=statename.value{3};
            %varStruct.(varname).string=['$' regexprep(statename.value{2},'([0-9]+)\>','_$1') '$'];
        end
        
        [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
        if verLessThan('symbolic','8')
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        else
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end
        end

    case 'EQUATIONNUM'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrievemodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
                varStruct.(varname).string{ii}=int2str(equationnum.value);
            end
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'EQUATIONNUMIMPLICIT'
        try
            statenum=retrievemodelinformation(ocStruct,'statenum');
            costatenum=retrievemodelinformation(ocStruct,'costatenum');
            totalalgebraicequationnum=retrievemodelinformation(ocStruct,'totalalgebraicequationnum');
            varStruct.(varname).string=int2str(statenum.value+costatenum.value+totalalgebraicequationnum.value);
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ALGEBRAICEQUATIONNUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            algebraicequationnum=0;
            counter=1;
            while counter<=numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{counter});
                equationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier{1});
                if equationnum.value
                    algebraicequationnum=1;
                    break
                end
                counter=counter+1;
            end
            varStruct.(varname).string=int2str(algebraicequationnum);
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPLICITVARIABLENUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationvariablenameimplicit=retrievemodelinformation(ocStruct,'equationvariablenameimplicit',arcidentifier);
                varStruct.(varname).string{ii}=int2str(length(equationvariablenameimplicit.value));
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'TOTALIMPLICITVARIABLENUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                totalalgebraicequationnum=retrievemodelinformation(ocStruct,'totalalgebraicequationnum');
                varStruct.(varname).string{ii}=int2str(totalalgebraicequationnum.value);
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPLICITNUM'
        % returns
        try
            implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
            varStruct.(varname).string=int2str(length(implicitnonlinearcontrolindex.value));
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMDIMENSION'
        % returns
        odeinfo=retrievemodelinformation(ocStruct,'odedim');
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcfield{ii});
                cansysdim=odeinfo.value+length(implicitnonlinearcontrolindex.value);
                varStruct.(varname).string{ii}=int2str(cansysdim);
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'IMPLICITINDEX'
        % returns
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                implicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcfield{ii});
                if ~isempty(implicitnonlinearcontrolindex.value)
                    if length(implicitnonlinearcontrolindex.value)==1
                        varStruct.(varname).string{ii}=int2str(implicitnonlinearcontrolindex.value);
                    else
                        idx=sprintf('%d,',implicitnonlinearcontrolindex.value);
                        idx(end)=[];
                        varStruct.(varname).string{ii}=['[' idx ']'];
                    end
                else
                    varStruct.(varname).string{ii}='[]';
                end
            end
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                %varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'ACTIVECONSTRAINTINDEX'
        % returns
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                nonzerolmmcindex=retrievemodelinformation(ocStruct,'nonzerolmmcindex',field2arcidentifier(arcfield{ii}));
                if ~isempty(nonzerolmmcindex.value)
                    if length(nonzerolmmcindex.value)==1
                        varStruct.(varname).string{ii}=int2str(nonzerolmmcindex.value);
                    else
                        idx=sprintf('%d,',nonzerolmmcindex.value);
                        idx(end)=[];
                        varStruct.(varname).string{ii}=['[' idx ']'];
                    end
                else
                    varStruct.(varname).string{ii}='[]';
                end
            end
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                %varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCALGEBRAICEQUATIONNUM'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier);
                varStruct.(varname).string{ii}=int2str(equationnum.value);
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCALGEBRAICEQUATIONNUMIMPLICIT'
        % returns
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrievemodelinformation(ocStruct,'totalalgebraicequationnum',arcidentifier);
                varStruct.(varname).string{ii}=int2str(equationnum.value);
            end
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'COSTATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=[int2str(statenum.value+1) ':' int2str(statenum.value+costatenum.value)];

    case 'IMPLICITCONTROLCOORD'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        for ii=1:arcnum.value
            fidx=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            if isempty(fidx.value)
                varStruct.(varname).string{ii}='[]';
            else
                if numel(fidx.value)>1
                    varStruct.(varname).string{ii}=['[' int2str(statenum.value+costatenum.value+(1:numel(fidx.value))) ']'];
                else
                    varStruct.(varname).string{ii}=int2str(statenum.value+costatenum.value+1);
                end
            end
        end

    case 'IMPLICITCONTROLINDEX'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            fidx=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            if isempty(fidx.value)
                varStruct.(varname).string{ii}='[]';
            else
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}=int2str(fidx.value);
                else
                    varStruct.(varname).string{ii}=['[' int2str(fidx.value) ']'];
                end
            end
        end

    case 'OPTIMALCONTROLDYNAMICSINDEX'
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        %totalimplicitnonlinearcontrolindex=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            fidx=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            idx=zeros(1,numel(fidx.value));
            for jj=1:numel(fidx.value)
                idx(jj)=jj;%find(totalimplicitnonlinearcontrolindex.value==fidx.value(jj));
            end
            if ~isempty(fidx.value)
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}=['DXDt(' int2str(statenum.value+costatenum.value+idx) ',:)=tmp(' int2str(fidx.value) ',:);'];
                else
                    varStruct.(varname).string{ii}=['DXDt([' int2str(statenum.value+costatenum.value+idx) '],:)=tmp([' int2str(fidx.value) '],:);'];
                end
            else
                varStruct.(varname).string{ii}='tmp=DXDt;';
            end
        end
        if arcnum.value>1
            varStruct.(varname).arcdependent=1;
        else
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string=varStruct.(varname).string{1};
        end
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;

    case 'OPTIMALCONTROLDYNAMICSJACOBIANINDEX'
        varStruct.(varname).arcdependent=1;
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        fidx=retrievemodelinformation(ocStruct,'totalimplicitvariableindex');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            %fidx=retrievemodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcidentifier.value{ii});
            idx=zeros(1,numel(fidx.value));
            for jj=1:numel(fidx.value)
                idx(jj)=jj;%find(totalimplicitnonlinearcontrolindex.value==fidx.value(jj));
            end
            if ~isempty(fidx.value)
                if length(fidx.value)==1
                    varStruct.(varname).string{ii}=['DDXDtDX(' int2str(statenum.value+costatenum.value+idx) ',:)=tmp(' int2str(fidx.value) ',:);'];
                else
                    varStruct.(varname).string{ii}=['DDXDtDX([' int2str(statenum.value+costatenum.value+idx) '],:)=tmp([' int2str(fidx.value) '],:);'];
                end
            else
                varStruct.(varname).string{ii}='tmp=DDXDtDX;';
            end
        end
        if arcnum.value>1
            varStruct.(varname).arcdependent=1;
        else
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string=varStruct.(varname).string{1};
        end
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;

    case 'ODEDIM'
        varStruct.(varname).arcdependent=1;
        odedim=retrievemodelinformation(ocStruct,'odedim');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        for ii=1:arcnum.value
            varStruct.(varname).string{ii}=int2str(odedim.value);
        end

    case 'AEDIM'
        varStruct.(varname).arcdependent=1;
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        for ii=1:arcnum.value
            algebraicequationnum=retrievemodelinformation(ocStruct,'algebraicequationnum',arcidentifier.value{ii});
            varStruct.(varname).string{ii}=int2str(algebraicequationnum.value);
        end

    case 'DAEORDER'
        if ~isfield(varStruct,'AEDIM')
            varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'AEDIM');
        end
        varStruct.(varname).arcdependent=1;
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        for ii=1:arcnum.value
            order=['[' int2str(ones(1,statenum.value+costatenum.value))];
            aedim=str2double(varStruct.AEDIM.string{ii});
            if aedim
                order=[order '  ' int2str(zeros(1, aedim))];
            end
            order=[order ']'];
            varStruct.(varname).string{ii}=order;
        end

    case 'PARAMETERVALUES'
        % the strings used for assigning parametervalues from the par
        % argument in functions to the spcific paraemter variable names,
        % e.g. r=pararg(1); ...
        varStruct.(varname).multline=2;
        basename=getbasicname('parametervariables');
        parvar=algebraicterm2string(ocStruct,'parametervariables');
        parvarexpr=retrievemodelinformation(ocStruct,'parametervalueexpression');
        for ii=1:numel(parvar.term)
            if isempty(parvarexpr.value{ii})
            varStruct.(varname).string{ii}=[parvar.term{ii} '=' basename '(' num2str(ii) ');'];
            else
                expr=parvarexpr.value{ii};
                for jj=1:numel(parvar.term)
                    if findstr(expr,parvar.term{jj})
                        expr=regexprep(expr,['\<' parvar.term{jj} '\>'],[basename '(' num2str(jj) ')']);
                    end
                end
                varStruct.(varname).string{ii}=[parvar.term{ii} '=' expr ';'];
            end
        end

    case 'ZEROSNTIMESNUMEQ'
        % returns a zero matrix where number of rows = number of states and
        % number of columns = number of equations of the canonical system
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            statenum=retrievemodelinformation(ocStruct,'statenum');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrievemodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
                zerovec=sprintf('%d,',zeros(1,equationnum.value));
                zerovec(end)=[];
                for jj=1:statenum.value
                    if jj==1
                        if statenum.value>1
                            varStruct.(varname).string{ii}{jj}=['[' zerovec '; ...'];
                        else
                            varStruct.(varname).string{ii}{jj}=['[' zerovec];
                        end
                    elseif jj<statenum.value
                        varStruct.(varname).string{ii}{jj}=[zerovec '; ...'];
                    else
                        varStruct.(varname).string{ii}{jj}=zerovec;
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} '];'];
            end
            varStruct.(varname).multline=2;
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ZEROSNUMEQ'
        % returns a zero matrix where number of rows = number of states and
        % number of columns = number of equations of the canonical system
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                equationnum=retrievemodelinformation(ocStruct,'canonicalsystemequationnum',arcidentifier{1});
                zerovec=sprintf('%d,',zeros(1,equationnum.value));
                zerovec(end)=[];
                varStruct.(varname).string{ii}=['[' zerovec ']'];
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).arcdependent=1;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEM'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'canonicalsystem',1,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end



    case 'NONSMOOTHFUNCTIONTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            for ii=1:numel(arcfield)
                nonsmoothfunctionterm=algebraicterm2string(ocStruct,'nonsmoothfunctionterm',1,arcfield{ii});
                for jj=1:numel(nonsmoothfunctionterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=nonsmoothfunctionterm.term{jj};
                    else
                        varStruct.(varname).string{ii}{jj}=nonsmoothfunctionterm.term{jj};
                    end
                end
                %varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'NONSMOOTHJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            nonsmoothfunctionnum=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                nonsmoothjacobianterm=algebraicterm2string(ocStruct,'nonsmoothjacobianterm',0,arcfield{ii});
                numlines=nonsmoothfunctionnum.value;
                counter=0;
                for jj=1:numel(nonsmoothjacobianterm.term)
                    if rem(jj,numlines)==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out(:,:,' num2str(counter) ')=' nonsmoothjacobianterm.term{jj}];
                    elseif rem(jj,numlines)==0
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out(:,:,' num2str(counter) ')=' nonsmoothjacobianterm.term{jj} ';'];
                    else
                        varStruct.(varname).string{ii}{jj}=nonsmoothjacobianterm.term{jj};
                    end
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        

    case 'NONSMOOTHPARAMETERJACOBIANTERM'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                nonsmoothparameterjacobianterm=algebraicterm2string(ocStruct,'nonsmoothparameterjacobianterm',0,arcfield{ii});
                for jj=1:numel(nonsmoothparameterjacobianterm.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' nonsmoothparameterjacobianterm.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=nonsmoothparameterjacobianterm.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'CANONICALSYSTEMDYNAMICS'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'canonicalsystem',1,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' canonicalsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMDYNAMICSCONT'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'canonicalsystemcont',1,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' canonicalsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMALGEBRAICEQUATION'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequation=algebraicterm2string(ocStruct,'algebraicequation',1,arcfield{ii});
                if ~isempty(algebraicequation.term)
                    for jj=1:numel(algebraicequation.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequation.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequation.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMALGEBRAICEQUATIONCONT'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequation=algebraicterm2string(ocStruct,'algebraicequationcont',1,arcfield{ii});
                if ~isempty(algebraicequation.term)
                    for jj=1:numel(algebraicequation.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequation.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequation.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ALGEBRAICEQUATIONIMPLICIT'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequationimplicit=algebraicterm2string(ocStruct,'algebraicequationimplicit',1,arcfield{ii});
                if ~isempty(algebraicequationimplicit.term)
                    for jj=1:numel(algebraicequationimplicit.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequationimplicit.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequationimplicit.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=3;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ALGEBRAICEQUATIONIMPLICITJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequationimplicitjacobian=algebraicterm2string(ocStruct,'algebraicequationimplicitjacobian',1,arcfield{ii});
                if ~isempty(algebraicequationimplicitjacobian.term)
                    for jj=1:numel(algebraicequationimplicitjacobian.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' algebraicequationimplicitjacobian.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequationimplicitjacobian.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=3;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ALGEBRAICEQUATIONIMPLICITPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                algebraicequationimplicitparameterjacobian=algebraicterm2string(ocStruct,'algebraicequationimplicitparameterjacobian',1,arcfield{ii});
                if ~isempty(algebraicequationimplicitparameterjacobian.term)
                    for jj=1:numel(algebraicequationimplicitparameterjacobian.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' algebraicequationimplicitparameterjacobian.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequationimplicitparameterjacobian.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=3;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'INVOPTIMALCONTROLDYNAMICSLEFTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'invoptimalcontroldynamicsleftside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsleftside.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLDYNAMICSLEFTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'optimalcontroldynamicsleftside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsleftside.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLDYNAMICS'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamics=algebraicterm2string(ocStruct,'optimalcontroldynamics',1,arcfield{ii});
                if ~isempty(optimalcontroldynamics.term)
                    for jj=1:numel(optimalcontroldynamics.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['dudt=' optimalcontroldynamics.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamics.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='dudt=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=3;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCONTROLDYNAMICS'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamics=algebraicterm2string(ocStruct,'symbolicoptimalcontroldynamics',0,arcfield{ii});
                if ~isempty(optimalcontroldynamics.term)
                    for jj=1:numel(optimalcontroldynamics.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['dudt=' optimalcontroldynamics.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamics.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='dudt=[];';
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'TENSOROPTIMALCONTROLDYNAMICSLEFTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'tensoroptimalcontroldynamicsleftside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    %varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLDYNAMICSRIGHTSIDE'
        try
            pontryaginfunctionDuDX=algebraicterm2string(ocStruct,'pontryaginfunctionDuDX',0);
            for jj=1:numel(pontryaginfunctionDuDX.term)
                if jj==1
                    varStruct.(varname).string{jj}=pontryaginfunctionDuDX.term{jj};
                else
                    varStruct.(varname).string{jj}=pontryaginfunctionDuDX.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONSTRAINEDOPTIMALCONTROLDYNAMICSRIGHTSIDE'
        arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
        arcnum=retrievemodelinformation(ocStruct,'arcnum');
        varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
        for ii=1:numel(arcfield)
            pontryaginfunctionDlmmcDX=algebraicterm2string(ocStruct,'pontryaginfunctionDlmmcDX',1,arcfield{ii});
            if ~isempty(pontryaginfunctionDlmmcDX.term)
                for jj=1:numel(pontryaginfunctionDlmmcDX.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['tmp=' pontryaginfunctionDlmmcDX.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=pontryaginfunctionDlmmcDX.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            else
                varStruct.(varname).string{ii}{1}='tmp=[];';
            end
        end
        varStruct.(varname).multline=1;
        if arcnum.value>1
            varStruct.(varname).arcdependent=1;
        else
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string=varStruct.(varname).string{1};
        end
        varStruct.(varname).type='algterm';

    case 'JACOBIANOPTIMALCONTROLDYNAMICSRIGHTSIDE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsleftside=algebraicterm2string(ocStruct,'jacobianoptimalcontroldynamicsrightside',0,arcfield{ii});
                if ~isempty(optimalcontroldynamicsleftside.term)
                    for jj=1:numel(optimalcontroldynamicsleftside.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsleftside.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsleftside.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='out=[];';
                end
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCANONICALSYSTEMDYNAMICS'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'symboliccanonicalsystem',0,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' canonicalsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICEQUILIBRIUMEQUATION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                equilibriumequation=algebraicterm2string(ocStruct,'equilibriumequation',0,arcfield{ii});
                for jj=1:numel(equilibriumequation.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' equilibriumequation.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=equilibriumequation.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCANONICALSYSTEMALGEBRAICEQUATION'
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            for ii=1:numel(arcfield)
                algebraicequation=algebraicterm2string(ocStruct,'symbolicalgebraicequation',0,arcfield{ii});
                if ~isempty(algebraicequation.term)
                    for jj=1:numel(algebraicequation.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['outae=' algebraicequation.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=algebraicequation.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                else
                    varStruct.(varname).string{ii}{1}='outae=mystr2sym([]);';
                end
            end
            varStruct.(varname).vectorize=0;
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'canonicalsystemjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECONTROLCANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'canonicalsystemjacobianstatecontrol',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATECOSTATEJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                statecostatejacobian=algebraicterm2string(ocStruct,'statecostatejacobian',0,arcfield{ii});
                for jj=1:numel(statecostatejacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' statecostatejacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=statecostatejacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SPECIFICCANONICALSYSTEMDYNAMICS'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                statecostatejacobian=algebraicterm2string(ocStruct,'specificcanonicalsystem',1,arcfield{ii});
                for jj=1:numel(statecostatejacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' statecostatejacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=statecostatejacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SPECIFICPONTRYAGINFUNCTION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                specificpontryaginfunction=algebraicterm2string(ocStruct,'specificpontryaginfunction',1,arcfield{ii});
                varStruct.(varname).string{ii}=['out=' specificpontryaginfunction.term ';'];
            end
            varStruct.(varname).multline=0;
            if 1%arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'symboliccanonicalsystemjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSTATECONTROLCANONICALSYSTEMJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'symbolicstatecontrolcanonicalsystemjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'SYMBOLICSTATECONTROLCANONICALSYSTEMDYNAMICS'
         try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemjacobian=algebraicterm2string(ocStruct,'symbolicstatecontrolcanonicalsystem',0,arcfield{ii});
                for jj=1:numel(canonicalsystemjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemparameterjacobian=algebraicterm2string(ocStruct,'canonicalsystemparameterjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemparameterjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemparameterjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemparameterjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'VARIATIONPARAMETERNUM'
        variationparameternum=retrievemodelinformation(ocStruct,'variationparameternum');
        varStruct.(varname).string=num2str(variationparameternum.value);


    case 'VARIATIONALSTATECOSTATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=[int2str(statenum.value+costatenum.value+1) ':' int2str(2*(statenum.value+costatenum.value))];

    case 'VARIATIONALSTATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=[int2str(statenum.value+costatenum.value+1) ':' int2str(2*statenum.value+costatenum.value)];

    case 'VARIATIONALCOSTATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatenum=retrievemodelinformation(ocStruct,'costatenum');
        varStruct.(varname).string=[int2str(2*statenum.value+costatenum.value+1) ':' int2str(2*(statenum.value+costatenum.value))];
        
    case 'VARIATIONALVARIABLE'
        try
            variationstatename=retrievemodelinformation(ocStruct,'variationstatename');
            variationstatenum=retrievemodelinformation(ocStruct,'variationstatenum');
            variationcostatename=retrievemodelinformation(ocStruct,'variationcostatename');
            variationcostatenum=retrievemodelinformation(ocStruct,'variationcostatenum');
            tmp=cellstr([char(variationstatename.value) char(repmat({';'},1,variationstatenum.value))]);
            tmp=[tmp{:}];
            tmp(isspace(tmp))=[];
            str=tmp;
            tmp=cellstr([char(variationcostatename.value) char(repmat({';'},1,variationcostatenum.value))]);
            tmp=[tmp{:}];
            tmp(isspace(tmp))=[];
            tmp(end)=[];
            str=['[' str tmp ']'];
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).string={str};
            varStruct.(varname).type='algterm';
            varStruct.(varname).multline=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end


    case 'VARIATIONALOPTIMALCONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'variationalcontrolvalue',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'VARIATIONALPARAMETEROBJECTIVEFUNCTION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                parameterderivativehamiltonianfunction=algebraicterm2string(ocStruct,'variationalparameterobjectivefunction',0,arcfield{ii});
                varStruct.(varname).string{ii}=['Jpar=' parameterderivativehamiltonianfunction.term{1} ';'];
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'VARIATIONALJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                variationaljacobian=algebraicterm2string(ocStruct,'variationaljacobian',0,arcfield{ii});
                for jj=1:numel(variationaljacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' variationaljacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=variationaljacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMVARIATIONPARAMETERDERIVATIVE'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemvariationparameterderivative=algebraicterm2string(ocStruct,'canonicalsystemvariationparameterderivative',0,arcfield{ii});
                for jj=1:numel(canonicalsystemvariationparameterderivative.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['Jpar=' canonicalsystemvariationparameterderivative.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemvariationparameterderivative.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'VARIATIONALHAMILTONIANFUNCTION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                variationalhamiltonianfunction=algebraicterm2string(ocStruct,'variationalhamiltonianfunction',0,arcfield{ii});
                varStruct.(varname).string{ii}=['HX=' variationalhamiltonianfunction.term{1} ';'];
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'VARIATIONALSALVAGEVALUE'
        try
            variationalsalvagevalue=algebraicterm2string(ocStruct,'variationalsalvagevalue',0);
            for jj=1:numel(variationalsalvagevalue.term)
                if jj==1
                    varStruct.(varname).string{jj}=['SX=' variationalsalvagevalue.term{jj}];
                else
                    varStruct.(varname).string{jj}=variationalsalvagevalue.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'PARAMETERVARIATIONALSALVAGEVALUE'
        try
            parametervariationalsalvagevalue=algebraicterm2string(ocStruct,'parametervariationalsalvagevalue',0);
            for jj=1:numel(parametervariationalsalvagevalue.term)
                if jj==1
                    varStruct.(varname).string{jj}=['SP=' parametervariationalsalvagevalue.term{jj}];
                else
                    varStruct.(varname).string{jj}=parametervariationalsalvagevalue.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONSTRAINTJACOBIAN'
        try
            constraintjacobian=algebraicterm2string(ocStruct,'constraintjacobian',0);
            for jj=1:numel(constraintjacobian.term)
                if jj==1
                    varStruct.(varname).string{jj}=['J=' constraintjacobian.term{jj}];
                elseif jj<numel(constraintjacobian.term)
                    varStruct.(varname).string{jj}=constraintjacobian.term{jj};
                else
                    varStruct.(varname).string{jj}=[constraintjacobian.term{jj}];
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONSTRAINTCONTROLJACOBIAN'
        try
            constraintjacobian=algebraicterm2string(ocStruct,'constraintcontroljacobian',0);
            for jj=1:numel(constraintjacobian.term)
                if jj==1
                    varStruct.(varname).string{jj}=['Ju=' constraintjacobian.term{jj}];
                elseif jj<numel(constraintjacobian.term)
                    varStruct.(varname).string{jj}=constraintjacobian.term{jj};
                else
                    varStruct.(varname).string{jj}=[constraintjacobian.term{jj}];
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONSTRAINTPARAMETERJACOBIAN'
        try
            constraintjacobian=algebraicterm2string(ocStruct,'constraintparameterjacobian',0);
            for jj=1:numel(constraintjacobian.term)
                if jj==1
                    varStruct.(varname).string{jj}=['Jpar=' constraintjacobian.term{jj}];
                elseif jj<numel(constraintjacobian.term)
                    varStruct.(varname).string{jj}=constraintjacobian.term{jj};
                else
                    varStruct.(varname).string{jj}=[constraintjacobian.term{jj}];
                end
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'PARAMETERDERIVATIVEHAMILTONIANFUNCTION'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                parameterderivativehamiltonianfunction=algebraicterm2string(ocStruct,'parameterderivativehamiltonianfunction',0,arcfield{ii});
                varStruct.(varname).string{ii}=['HP=' parameterderivativehamiltonianfunction.term{1} ';'];
            end
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

 
    case 'VARIATIONALTRANSVERSALITYCONDITION'
        try
            variationaltransversalitycondition=algebraicterm2string(ocStruct,'variationaltransversalitycondition',0);
            for jj=1:numel(variationaltransversalitycondition.term)
                if jj==1
                    varStruct.(varname).string{jj}=['SxX=' variationaltransversalitycondition.term{jj}];
                else
                    varStruct.(varname).string{jj}=variationaltransversalitycondition.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
 
    case 'PARAMETERDERIVATIVETRANSVERSALITYCONDITION'
        try
            parameterderivativetransversalitycondition=algebraicterm2string(ocStruct,'parameterderivativetransversalitycondition',0);
            for jj=1:numel(parameterderivativetransversalitycondition.term)
                if jj==1
                    varStruct.(varname).string{jj}=['SxP=' parameterderivativetransversalitycondition.term{jj}];
                else
                    varStruct.(varname).string{jj}=parameterderivativetransversalitycondition.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'VARIATIONALSTATECONSTRAINT'
        try
            inequalitystateconstraint=algebraicterm2string(ocStruct,'variationalstateconstraint',1,[],'inequality');
            for jj=1:numel(inequalitystateconstraint.term)
                if jj==1
                    varStruct.(varname).string{jj}=['SCx=[' inequalitystateconstraint.term{jj}];
                    if numel(inequalitystateconstraint.term)>1
                        varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '; ...'];
                    end
                elseif jj<numel(inequalitystateconstraint.term)
                    varStruct.(varname).string{jj}=[inequalitystateconstraint.term{jj} '; ...'];
                else
                    varStruct.(varname).string{jj}=inequalitystateconstraint.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '];'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
       
    case 'STATECONTROLCANONICALSYSTEMPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemparameterjacobian=algebraicterm2string(ocStruct,'statecontrolcanonicalsystemparameterjacobian',0,arcfield{ii});
                for jj=1:numel(canonicalsystemparameterjacobian.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemparameterjacobian.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemparameterjacobian.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMDERIVATIVETIME'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemderivativetime=algebraicterm2string(ocStruct,'canonicalsystemderivativetime',0,arcfield{ii});
                for jj=1:numel(canonicalsystemderivativetime.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' canonicalsystemderivativetime.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemderivativetime.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCIDENTIFIERALGEBRAICEQUATION'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'algebraicequation');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'ARCIDENTIFIERCANONICALSYSTEMJACOBIAN'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'CANONICALSYSTEMHESSIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemhessian=algebraicterm2string(ocStruct,'canonicalsystemhessian',0,arcfield{ii});
                numlines=canonicalsystemhessian.info;
                counter=0;
                for jj=1:numel(canonicalsystemhessian.term)
                    if rem(jj,numlines)==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out(:,:,' num2str(counter) ')=' canonicalsystemhessian.term{jj}];
                    elseif rem(jj,numlines)==0
                        varStruct.(varname).string{ii}{jj}=[canonicalsystemhessian.term{jj} ';'];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemhessian.term{jj};
                    end
                end
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMTOTALHESSIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                canonicalsystemparameterhessian=algebraicterm2string(ocStruct,'canonicalsystemtotalhessian',0,arcfield{ii});
                numlines=canonicalsystemparameterhessian.info;
                counter=0;
                for jj=1:numel(canonicalsystemparameterhessian.term)
                    if rem(jj,numlines)==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out(:,:,' num2str(counter) ')=' canonicalsystemparameterhessian.term{jj}];
                    elseif rem(jj,numlines)==0
                        varStruct.(varname).string{ii}{jj}=[canonicalsystemparameterhessian.term{jj} ';'];
                    else
                        varStruct.(varname).string{ii}{jj}=canonicalsystemparameterhessian.term{jj};
                    end
                end
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2PONTRYAGINFUNCTIONDU2'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcontroldynamicsPuu=algebraicterm2string(ocStruct,'optimalcontroldynamicsPuu',1,arcfield{ii});
                counter=0;
                for jj=1:numel(optimalcontroldynamicsPuu.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' optimalcontroldynamicsPuu.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=optimalcontroldynamicsPuu.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPLICITGX'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                implicitcontroldynamicsGX=algebraicterm2string(ocStruct,'implicitcontroldynamicsGX',0,arcfield{ii});
                counter=0;
                for jj=1:numel(implicitcontroldynamicsGX.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' implicitcontroldynamicsGX.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=implicitcontroldynamicsGX.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'IMPLICITG'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                implicitcontroldynamicsfunction=algebraicterm2string(ocStruct,'implicitcontroldynamicsfunction',0,arcfield{ii});
                counter=0;
                for jj=1:numel(implicitcontroldynamicsfunction.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' implicitcontroldynamicsfunction.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=implicitcontroldynamicsfunction.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'IMPLICITGP'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                implicitcontroldynamicsGP=algebraicterm2string(ocStruct,'implicitcontroldynamicsGP',0,arcfield{ii});
                counter=0;
                for jj=1:numel(implicitcontroldynamicsGP.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' implicitcontroldynamicsGP.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=implicitcontroldynamicsGP.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'IMPLICITGU'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                implicitcontroldynamicsGU=algebraicterm2string(ocStruct,'implicitcontroldynamicsGU',0,arcfield{ii});
                counter=0;
                for jj=1:numel(implicitcontroldynamicsGU.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' implicitcontroldynamicsGU.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=implicitcontroldynamicsGU.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DLAGRANGIANDU'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                dlagrangiandU=algebraicterm2string(ocStruct,'dlagrangiandU',1,arcfield{ii});
                counter=0;
                for jj=1:numel(dlagrangiandU.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' dlagrangiandU.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=dlagrangiandU.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2LAGRANGIANDUI2'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                d2lagrangiandui2=algebraicterm2string(ocStruct,'d2lagrangiandui2',1,arcfield{ii});
                counter=0;
                for jj=1:numel(d2lagrangiandui2.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' d2lagrangiandui2.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=d2lagrangiandui2.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2LAGRANGIANDU2'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                d2lagrangiandU2=algebraicterm2string(ocStruct,'d2lagrangiandU2',1,arcfield{ii});
                counter=0;
                for jj=1:numel(d2lagrangiandU2.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' d2lagrangiandU2.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=d2lagrangiandU2.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2LAGRANGIANDUDX'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                d2lagrangiandUdX=algebraicterm2string(ocStruct,'d2lagrangiandUdX',1,arcfield{ii});
                counter=0;
                for jj=1:numel(d2lagrangiandUdX.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' d2lagrangiandUdX.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=d2lagrangiandUdX.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2LAGRANGIANDUDP'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                d2lagrangiandUdP=algebraicterm2string(ocStruct,'d2lagrangiandUdP',0,arcfield{ii});
                counter=0;
                for jj=1:numel(d2lagrangiandUdP.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' d2lagrangiandUdP.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=d2lagrangiandUdP.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2HAMILTONIANDU2'
        try
            pontryaginfunctionDu2=algebraicterm2string(ocStruct,'pontryaginfunctionDu2');
            counter=0;
            for jj=1:numel(pontryaginfunctionDu2.term)
                if jj==1
                    counter=counter+1;
                    varStruct.(varname).string{jj}=['out=' pontryaginfunctionDu2.term{jj}];
                else
                    varStruct.(varname).string{jj}=pontryaginfunctionDu2.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2HAMILTONIANDX2'
        try
            D2pontryaginfunctionDX2=algebraicterm2string(ocStruct,'D2pontryaginfunctionDX2');
            counter=0;
            for jj=1:numel(D2pontryaginfunctionDX2.term)
                if jj==1
                    counter=counter+1;
                    varStruct.(varname).string{jj}=['out=' D2pontryaginfunctionDX2.term{jj}];
                else
                    varStruct.(varname).string{jj}=D2pontryaginfunctionDX2.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2HAMILTONIANDUI2'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                D2pontryaginfunctionDui2=algebraicterm2string(ocStruct,'D2pontryaginfunctionDui2',0,arcfield{ii});
                counter=0;
                for jj=1:numel(D2pontryaginfunctionDui2.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' D2pontryaginfunctionDui2.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=D2pontryaginfunctionDui2.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2HAMILTONIANDUIDX'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                D2pontryaginfunctionDUDX=algebraicterm2string(ocStruct,'D2pontryaginfunctionDUDX',0,arcfield{ii});
                counter=0;
                for jj=1:numel(D2pontryaginfunctionDUDX.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' D2pontryaginfunctionDUDX.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=D2pontryaginfunctionDUDX.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'D2HAMILTONIANDUDX'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                D2pontryaginfunctionDUDX=algebraicterm2string(ocStruct,'D2pontryaginfunctionDUDX',0,arcfield{ii});
                counter=0;
                for jj=1:numel(D2pontryaginfunctionDUDX.term)
                    if jj==1
                        counter=counter+1;
                        varStruct.(varname).string{ii}{jj}=['out=' D2pontryaginfunctionDUDX.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=D2pontryaginfunctionDUDX.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'PARAMETERHESSIANCOORD'
        arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
        varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
        parameternum=retrievemodelinformation(ocStruct,'parameternum');
        for ii=1:numel(arcfield)
            equationnum=retrievemodelinformation(ocStruct,'canonicalsystemequationnum',field2arcidentifier(arcfield{ii}));
            varStruct.(varname).string{ii}=['1:' int2str(equationnum.value) ',1:'  int2str(equationnum.value) ',' int2str(equationnum.value+1) ':' int2str(parameternum.value+equationnum.value)];
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        end

    case 'ARCIDENTIFIERLAGRANGEMULTIPLIER'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'lagrangemultiplier');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'OPTIMALCONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'controlvalue',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLVALUEIMPLICIT'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'controlvalue',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCONTROLVALUE4STATECONTROL'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'optimalvalue4statecontrol',1,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OPTIMALCOSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                optimalcostatevalue=algebraicterm2string(ocStruct,'optimalcostatevalue',1,arcfield{ii});
                for jj=1:numel(optimalcostatevalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' optimalcostatevalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=optimalcostatevalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'EXPLICITSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                explicitstatevalue=algebraicterm2string(ocStruct,'explicitstatevalue',1,arcfield{ii});
                for jj=1:numel(explicitstatevalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' explicitstatevalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=explicitstatevalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                symbolicstatevalue=algebraicterm2string(ocStruct,'symbolicstatevalue',0,arcfield{ii});
                for jj=1:numel(symbolicstatevalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' symbolicstatevalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=symbolicstatevalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCONTROLVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'symboliccontrolvalue',0,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCONTROLVALUE4STATECONTROL'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'symboliccontrolvalue4statecontrol',0,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOPTIMALCOSTATEVALUE'
        try
            arcfield=getarcclass(ocStruct,'controlvalue');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                controlvalue=algebraicterm2string(ocStruct,'symboliccostatevalue',0,arcfield{ii});
                for jj=1:numel(controlvalue.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' controlvalue.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=controlvalue.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONTROLLAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultcc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            for ii=1:numel(arcid)
                lagrangemultcc=algebraicterm2string(ocStruct,'lagrangemultcc',1,arcid{ii});
                for jj=1:numel(lagrangemultcc.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' lagrangemultcc.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=lagrangemultcc.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATELAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultsc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            if isempty(arcid)
                varStruct.(varname).string='';
            else
                for ii=1:numel(arcid)
                    lagrangemultsc=algebraicterm2string(ocStruct,'lagrangemultsc',1,arcid{ii});
                    for jj=1:numel(lagrangemultsc.term)
                        if jj==1
                            varStruct.(varname).string{ii}{jj}=['out=' lagrangemultsc.term{jj}];
                        else
                            varStruct.(varname).string{ii}{jj}=lagrangemultsc.term{jj};
                        end
                    end
                    varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
                end
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %%lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICCONTROLLAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultcc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            for ii=1:numel(arcid)
                lagrangemultcc=algebraicterm2string(ocStruct,'symboliclagrangemultcc',0,arcid{ii});
                for jj=1:numel(lagrangemultcc.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' lagrangemultcc.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=lagrangemultcc.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICSTATELAGRANGEMULTIPLIER'
        try
            arcid=getarcclass(ocStruct,'lagrangemultsc');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcid);
            for ii=1:numel(arcid)
                lagrangemultsc=algebraicterm2string(ocStruct,'symboliclagrangemultsc',0,arcid{ii});
                for jj=1:numel(lagrangemultsc.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' lagrangemultsc.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=lagrangemultsc.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).arcdependent=1;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'ARCIDENTIFIERCONTROL'
        varStruct.(varname).arcdependent=1;
        try
            arcfield=getarcclass(ocStruct,'control');
            for ii=1:numel(arcfield)
                arcidentifier=field2arcidentifier(arcfield{ii});
                if numel(arcidentifier)>1
                    varStruct.(varname).string{ii}=['{' arcidentifier{1}];
                    for jj=2:numel(arcidentifier)
                        varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ',' arcidentifier{jj}];
                    end
                    varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} '}'];
                else
                    varStruct.(varname).string{ii}=arcidentifier{1};
                end
            end
        catch
            %lasterr
            varStruct.(varname).string={'0'};
        end

    case 'SYMBOLICDISCOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'symbolicdiscobjectivefunction',0);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'discobjectivefunction',1);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTIONDERIVATIVETIME'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:numel(arcfield)
                discobjectivefunctionderivativetime=algebraicterm2string(ocStruct,'discobjectivefunctionderivativetime',0,arcfield{ii});
                for jj=1:numel(discobjectivefunctionderivativetime.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' discobjectivefunctionderivativetime.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=discobjectivefunctionderivativetime.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTIONJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                discobjectivefunctionDx=algebraicterm2string(ocStruct,'discobjectivefunctionjacobian',0,arcfield{ii});
                for jj=1:numel(discobjectivefunctionDx.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' discobjectivefunctionDx.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=discobjectivefunctionDx.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOBJECTIVEFUNCTIONPARAMETERJACOBIAN'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                discobjectivefunctionDx=algebraicterm2string(ocStruct,'discobjectivefunctionparameterjacobian',0,arcfield{ii});
                for jj=1:numel(discobjectivefunctionDx.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' discobjectivefunctionDx.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=discobjectivefunctionDx.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'HAMILTONIANFUNCTION'
        try
            hamiltonianfunction=algebraicterm2string(ocStruct,'hamiltonianfunction',1);
            varStruct.(varname).string{1}=['out=' hamiltonianfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'HAMILTONIAN4GRAD'
        try
            hamiltonianfunction=algebraicterm2string(ocStruct,'hamiltonianfunction',1);
            varStruct.(varname).string{1}=['out=' hamiltonianfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'PONTRYAGINFUNCTION'
        try
            pontryaginfunction=algebraicterm2string(ocStruct,'pontryaginfunction',1);
            varStruct.(varname).string{1}=['out=' pontryaginfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SALVAGEVALUE'
        try
            salvagevalue=algebraicterm2string(ocStruct,'salvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DISCOUNTRATE'
        try
            discountratevariable=retrievemodelinformation(ocStruct,'discountratevariable');
            varStruct.(varname).string{1}=discountratevariable.value;
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'DISCOUNTEDSALVAGEVALUE'
        try
            salvagevalue=algebraicterm2string(ocStruct,'discountedsalvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICOBJECTIVEFUNCTION'
        try
            symbolicobjectivefunction=algebraicterm2string(ocStruct,'symbolicobjectivefunction',0);
            varStruct.(varname).string{1}=['out=' symbolicobjectivefunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'SYMBOLICPONTRYAGINFUNCTION'
        try
            symbolicpontryaginfunction=algebraicterm2string(ocStruct,'symbolicpontryaginfunction',0);
            varStruct.(varname).string{1}=['out=' symbolicpontryaginfunction.term{1} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DHAMILTONIANDX'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                pontryaginfunctionDx=algebraicterm2string(ocStruct,'pontryaginfunctionDX',0,arcfield{ii});
                for jj=1:numel(pontryaginfunctionDx.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' pontryaginfunctionDx.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=pontryaginfunctionDx.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DHAMILTONIANDU'
        try
            arcfield=getarcclass(ocStruct,'canonicalsystemjacobian');
            arcnum=retrievemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).arcidentifier=field2arcidentifier(arcfield);
            for ii=1:arcnum.value
                %DpontryaginfunctionDu=algebraicterm2string(ocStruct,'DpontryaginfunctionDu',0,arcfield{ii});
                DpontryaginfunctionDu=algebraicterm2string(ocStruct,'DHamiltonianDu',1,arcfield{ii});
                for jj=1:numel(DpontryaginfunctionDu.term)
                    if jj==1
                        varStruct.(varname).string{ii}{jj}=['out=' DpontryaginfunctionDu.term{jj}];
                    else
                        varStruct.(varname).string{ii}{jj}=DpontryaginfunctionDu.term{jj};
                    end
                end
                varStruct.(varname).string{ii}{jj}=[varStruct.(varname).string{ii}{jj} ';'];
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            if arcnum.value>1
                varStruct.(varname).arcdependent=1;
            else
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).string=varStruct.(varname).string{1};
            end
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'STATEDYNAMICS'
        try
            statedynamics=algebraicterm2string(ocStruct,'statedynamics',1,'');
            for jj=1:numel(statedynamics.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' statedynamics.term{jj}];
                else
                    varStruct.(varname).string{jj}=statedynamics.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='bvp';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

 
    case 'ADJOINTDYNAMICS'
        try
            adjointsystem=algebraicterm2string(ocStruct,'adjointdynamics',1,'');
            for jj=1:numel(adjointsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' adjointsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=adjointsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'TRANSVERSALITYCONDITION'
        try
            transversalitycondition=algebraicterm2string(ocStruct,'transversalitycondition',0);
            for jj=1:numel(transversalitycondition.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' transversalitycondition.term{jj}];
                else
                    varStruct.(varname).string{jj}=transversalitycondition.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='bvp';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'OBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'objectivefunction',1);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DATAPATHNAME'
        varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct))) ''''];

    case 'RESULTFILE'
        varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct)),'SaveIntermediateResults') ''''];

    case 'GLOBALVARFILE'
        varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct)),'SaveIntermediateResultsGlobalVariable') ''''];

    case 'INFODETAILS'
        varStruct.(varname).multline=2;
        year=datevec(now);
        year=year(1);
        varStruct.(varname).string={'%', ...
            ['% this file was automatically created: ' datestr(now)], ...
            ['% written by Dieter Grass, 2001 - ' num2str(year)]};
        
%%%%%%
%% for the gradient method        
%%%%%

    case 'LOCSTATE'
        varStruct.(varname).string=getbasicname('locstate');
        varStruct.(varname).method='grad';

    case 'LOCCOSTATE'
        varStruct.(varname).string=getbasicname('loccostate');
        varStruct.(varname).method='grad';

    case 'LOCCONTROL'
        varStruct.(varname).string=getbasicname('loccontrol');
        varStruct.(varname).method='grad';

    case 'LOCSTATECONTROL'
        varStruct.(varname).string='loc_statecontrol';
        varStruct.(varname).method='grad';

    case 'GRADIENTLOCCONTROL'
        varStruct.(varname).string=getbasicname('gradientloccontrol');
        varStruct.(varname).method='grad';

    case 'LOCSTATECOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        if statenum.value==1
            varStruct.(varname).string=int2str(statenum.value);
        else
            varStruct.(varname).string=['1:' int2str(statenum.value)];
        end
        varStruct.(varname).method='grad';

    case 'LOCCONTROLNUM'
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        varStruct.(varname).string=int2str(controlnum.value);
        varStruct.(varname).method='grad';


    case 'LOCCONTROLCOORD'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        if controlnum.value==1
            varStruct.(varname).string=int2str(statenum.value+1);
        else
            varStruct.(varname).string=[int2str(statenum.value+1) ':' int2str(statenum.value+controlnum.value)];
        end
        varStruct.(varname).method='grad';

    case 'PARGAMMA'
        varStruct.(varname).string='par_gamma';
        varStruct.(varname).method='grad';

    case 'LOCALSTATEDYNAMICS'
        try
            statedynamics=algebraicterm2string(ocStruct,'localstatedynamics',1,'');
            for jj=1:numel(statedynamics.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' statedynamics.term{jj}];
                else
                    varStruct.(varname).string{jj}=statedynamics.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'LOCALCOSTATEDYNAMICS'
        try
            adjointsystem=algebraicterm2string(ocStruct,'localcostatedynamics',1,'');
            for jj=1:numel(adjointsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' adjointsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=adjointsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'LOCALOBJECTIVEFUNCTION'
        try
            objectivefunction=algebraicterm2string(ocStruct,'discobjectivefunction',1);
            varStruct.(varname).string{1}=['out=' objectivefunction.term{1} ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DHAMILTONIANDU'
        try
            hamiltonianDu=algebraicterm2string(ocStruct,'hamiltonianDu',1);
            for ii=1:numel(hamiltonianDu.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' hamiltonianDu.term{ii}];
                else
                    varStruct.(varname).string{ii}=hamiltonianDu.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'LOCALSALVAGEVALUE'
        try
            salvagevalue=algebraicterm2string(ocStruct,'discountedsalvagevalue',1);
            varStruct.(varname).string{1}=['out=' salvagevalue.term ';'];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'LOCALTRANSVERSALITYCONDITION'
        try
            transversalitycondition=algebraicterm2string(ocStruct,'transversalitycondition',0);
            for jj=1:numel(transversalitycondition.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' transversalitycondition.term{jj}];
                else
                    varStruct.(varname).string{jj}=transversalitycondition.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'KUHNTUCKERPOINT'
        try
            kuhntuckersolution=algebraicterm2string(ocStruct,'kuhntuckersolution',1);
            varStruct.(varname).string=kuhntuckersolution.term;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
        
    case 'KUHNTUCKERNUM'
        try
            kuhntuckersolution=retrievemodelinformation(ocStruct,'kuhntuckersolution','',getsymkernel());
            varStruct.(varname).string=int2str(length(kuhntuckersolution.value));
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
    case 'LOCALCONTROLPROJECTION'
        try
            locprojection=algebraicterm2string(ocStruct,'locprojection',1);
            varStruct.(varname).string=locprojection.term;
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
    case 'LOCALCONSTRAINT'
        try
            constraint=algebraicterm2string(ocStruct,'constraint',1);
            varStruct.(varname).string='';
            for jj=1:numel(constraint.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=[' constraint.term{jj}];
                    if numel(constraint.term)>1
                        varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '; ...'];
                    end
                elseif jj<numel(constraint.term)
                    varStruct.(varname).string{jj}=[constraint.term{jj} '; ...'];
                else
                    varStruct.(varname).string{jj}=constraint.term{jj};
                end
            end
            if ~isempty(varStruct.(varname).string)
                varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} '];'];
            else
                varStruct.(varname).string{1}='out=[];';
            end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='grad';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        
        %DAE
        
    case 'CANONICALSYSTEMDYNAMICSDAE'
        try
            canonicalsystem=algebraicterm2string(ocStruct,'canonicalsystem',1,'');
            for jj=1:numel(canonicalsystem.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' canonicalsystem.term{jj}];
                else
                    varStruct.(varname).string{jj}=canonicalsystem.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            %end
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end
        varStruct.(varname).method='dae';

    case 'CANONICALSYSTEMJACOBIANDAE'
        try
            canonicalsystemjacobiandae=algebraicterm2string(ocStruct,'canonicalsystemjacobiandae',0);
            for ii=1:numel(canonicalsystemjacobiandae.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' canonicalsystemjacobiandae.term{ii}];
                else
                    varStruct.(varname).string{ii}=canonicalsystemjacobiandae.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CANONICALSYSTEMPARAMETERJACOBIANDAE'
        try
            canonicalsystemparameterjacobiandae=algebraicterm2string(ocStruct,'canonicalsystemparameterjacobiandae',0);
            for ii=1:numel(canonicalsystemparameterjacobiandae.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' canonicalsystemparameterjacobiandae.term{ii}];
                else
                    varStruct.(varname).string{ii}=canonicalsystemparameterjacobiandae.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DLAGRANGEFUNCTIONDU'
        try
            dlagrangefunctiondu=algebraicterm2string(ocStruct,'dlagrangefunctiondu',1);
            for ii=1:numel(dlagrangefunctiondu.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' dlagrangefunctiondu.term{ii}];
                else
                    varStruct.(varname).string{ii}=dlagrangefunctiondu.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DLAGRANGEFUNCTIONDUJACOBIAN'
        try
            dlagrangefunctiondujacobian=algebraicterm2string(ocStruct,'dlagrangefunctiondujacobian',0);
            for ii=1:numel(dlagrangefunctiondujacobian.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' dlagrangefunctiondujacobian.term{ii}];
                else
                    varStruct.(varname).string{ii}=dlagrangefunctiondujacobian.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'DLAGRANGEFUNCTIONDUPARAMETERJACOBIAN'
        try
            dlagrangefunctionduparameterjacobian=algebraicterm2string(ocStruct,'dlagrangefunctionduparameterjacobian',0);
            for ii=1:numel(dlagrangefunctionduparameterjacobian.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' dlagrangefunctionduparameterjacobian.term{ii}];
                else
                    varStruct.(varname).string{ii}=dlagrangefunctionduparameterjacobian.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OBJECTIVEFUNCTIONJACOBIANDAE'
        try
            objectivefunctionjacobiandae=algebraicterm2string(ocStruct,'objectivefunctionjacobiandae',1);
            varStruct.(varname).string{1}=['out=' objectivefunctionjacobiandae.term{1}];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'OBJECTIVEFUNCTIONPARAMETERJACOBIANDAE'
        try
            objectivefunctionparameterjacobiandae=algebraicterm2string(ocStruct,'objectivefunctionparameterjacobiandae',1);
            varStruct.(varname).string{1}=['out=' objectivefunctionparameterjacobiandae.term{1}];
            varStruct.(varname).multline=0;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).type='algterm';
            varStruct.(varname).arcdependent=0;
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONTROLCOORDINATE'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        if controlnum.value==1
            varStruct.(varname).string=int2str(2*statenum.value+1);
        else
            varStruct.(varname).string=[int2str(2*statenum.value+1) ':' int2str(2*statenum.value+controlnum.value)];
        end
        varStruct.(varname).method='dae';

    case 'LAGRANGEMULTCCCOORDINATE'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if inequalitycontrolconstraintnum.value==1
            varStruct.(varname).string=int2str(2*statenum.value+controlnum.value+1);
        else
            varStruct.(varname).string=[int2str(2*statenum.value+controlnum.value+1) ':' int2str(2*statenum.value+controlnum.value+inequalitycontrolconstraintnum.value)];
        end
        varStruct.(varname).method='dae';

    case 'MIXEDCONTROLCONSTRAINT'
        mixedinequalityconstraintnum=retrievemodelinformation(ocStruct,'mixedinequalityconstraint','',getsymkernel());
        varStruct.(varname).string=num2str(mixedinequalityconstraintnum.value);
        varStruct.(varname).method='dae';
        
    case 'LAGRANGEMULTCCJACOBIANDAE'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        
        varStruct.(varname).string=['Jm=[zeros(' int2str(inequalitycontrolconstraintnum.value) ',' int2str(2*statenum.value+controlnum.value) ') eye(' int2str(inequalitycontrolconstraintnum.value) ')];'];
        varStruct.(varname).method='dae';

    case 'CONTROLCONSTRAINTJACOBIANDAE'
        try
            controlconstraintjacobiandae=algebraicterm2string(ocStruct,'controlconstraintjacobiandae');
            for ii=1:numel(controlconstraintjacobiandae.term)
                if ii==1
                    varStruct.(varname).string{ii}=['Jc=' controlconstraintjacobiandae.term{ii}];
                else
                    varStruct.(varname).string{ii}=controlconstraintjacobiandae.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'CONTROLCONSTRAINTPARAMETERJACOBIANDAE'
        try
            controlconstraintparameterjacobiandae=algebraicterm2string(ocStruct,'controlconstraintparameterjacobiandae');
            for ii=1:numel(controlconstraintparameterjacobiandae.term)
                if ii==1
                    varStruct.(varname).string{ii}=['Jc=' controlconstraintparameterjacobiandae.term{ii}];
                else
                    varStruct.(varname).string{ii}=controlconstraintparameterjacobiandae.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'TRANSVERSALITYCONDITIONDAE'
        try
            transversalityconditiondae=algebraicterm2string(ocStruct,'transversalityconditiondae',0);
            for jj=1:numel(transversalityconditiondae.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' transversalityconditiondae.term{jj}];
                else
                    varStruct.(varname).string{jj}=transversalityconditiondae.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).type='algterm';
            varStruct.(varname).method='dae';
        catch
            %lasterr
            varStruct.(varname).string=varname;
        end

    case 'BURMEISTERFISCHERMU'
        varStruct.(varname).string='bfmu';

    case 'BURMEISTERFISCHEREPSILON'
        varStruct.(varname).string='bfepsilon';
        
    case 'STATENAME'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        statename=retrievemodelinformation(ocStruct,'statename');
        tmp='';
        for ii=1:statenum.value
            if ii<statenum.value
                tmp=[tmp '$' statename.value{ii} '$,\ '];
            else
                tmp=[tmp '$' statename.value{ii} '$'];
            end
        end
        tmp=[tmp ''];
        varStruct.(varname).string{1}=tmp;
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
    case 'COSTATENAME'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        costatename=retrievemodelinformation(ocStruct,'costatename');
        tmp='';
        for ii=1:statenum.value
            if ii<statenum.value
                tmp=[tmp '$' costatename.value{ii} '$,\ '];
            else
                tmp=[tmp '$' costatename.value{ii} '$'];
            end
        end
        tmp=[tmp ''];
        varStruct.(varname).string{1}=tmp;
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
    case 'CONTROLNAME'
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        controlname=retrievemodelinformation(ocStruct,'controlname');
        tmp='';
        for ii=1:controlnum.value
            if ii<controlnum.value
                tmp=[tmp '$' controlname.value{ii} '$,\ '];
            else
                tmp=[tmp '$' controlname.value{ii} '$'];
            end
        end
        tmp=[tmp ''];
        varStruct.(varname).string{1}=tmp;
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
    case 'LAGRANGEMULTCCNAME'
        lagrangemultipliercontrolname=retrievemodelinformation(ocStruct,'lagrangemultipliercontrolname');
        tmp='';
        for ii=1:length(lagrangemultipliercontrolname.value)
            if ii<length(lagrangemultipliercontrolname.value)
                tmp=[tmp '$' lagrangemultipliercontrolname.value{ii} '$,\ '];
            else
                tmp=[tmp '$' lagrangemultipliercontrolname.value{ii} '$'];
            end
        end
        tmp=[tmp ''];
        varStruct.(varname).string{1}=tmp;
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
    case 'INDEPENDENTNAME'
        varStruct.(varname).string{1}=[ '$' getbasicname('independent') '$'];
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
        
    case 'ONESCOMPONENTNUMBER'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        varStruct.(varname).string{1}=[ '[' int2str(ones(1,2*statenum.value+inequalitycontrolconstraintnum.value+controlnum.value)) ']'];
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
    case 'TWOSCOMPONENTNUMBER'
        statenum=retrievemodelinformation(ocStruct,'statenum');
        controlnum=retrievemodelinformation(ocStruct,'controlnum');
        inequalitycontrolconstraintnum=retrievemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        varStruct.(varname).string{1}=[ '[' int2str(2*ones(1,2*statenum.value+inequalitycontrolconstraintnum.value+controlnum.value)) ']'];
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        
    case 'ONESPARAMETERNUMBER'
        parameternum=retrievemodelinformation(ocStruct,'parameternum');
        varStruct.(varname).string{1}=[ '[' int2str(ones(1,parameternum.value)) ']'];
        varStruct.(varname).multline=0;
        varStruct.(varname).vectorize=0;
        varStruct.(varname).type='algterm';
        varStruct.(varname).method='dae';
        

end