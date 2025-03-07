function varStruct=generatevariablestruct4odemodel(ocStruct,varStruct,varname,overwrite,varargin)
%
% a structure is generated which is then used to replace variables in
% template model files
% with overwrite 1 an already existing field of varStruct is overwritten
% the structure consists of a field given by the name of the variable. If
% the actual value of the variable depends on the arc the field string
% contains a cell array of strings and the field arcependent is set to one.
% multline 1 ... lines>1 are intended, e.g., DYNAMICSDYNAMICS
% multline 2 ...lines are not intended, e.g., PARAMETERVALUES
if nargin<=3
    overwrite=0;
end
if ~overwrite
    if isfield(varStruct,varname)
        return
    end
end
varStruct.(varname).arcdependent=0;
varStruct.(varname).arcidentifier=[];
varStruct.(varname).vectorize=0;
varStruct.(varname).type='';
varStruct.(varname).multline=0;

try
    switch varname
        case 'INDEPENDENT'
            if isfield(ocStruct.variable.independent,'name') && ~isempty(ocStruct.variable.independent.name)
                varStruct.(varname).string=ocStruct.variable.independent.name;
            else
                varStruct.(varname).string=getbasicname('independent');
            end
            
        case 'DEPENDENTVAR'
            varStruct.(varname).string=getbasicname('dependent');

        case 'PARVAR'
            varStruct.(varname).string=getbasicname('parametervariables');

        case 'ARCVAR'
            varStruct.(varname).string=getbasicname('arcidentifiervar');

        case 'ARCIDENTIFIER'
            varStruct.(varname).arcdependent=1;
            arcidentifier=retrieveodemodelinformation(ocStruct,'arcidentifier');
            varStruct.(varname).string=arcidentifier.value;

        case 'ODEDIM'
            varStruct.(varname).arcdependent=1;
            odedim=retrieveodemodelinformation(ocStruct,'odedim');
            arcnum=retrieveodemodelinformation(ocStruct,'arcnum');
            for ii=1:arcnum.value
                varStruct.(varname).string{ii}=int2str(odedim.value);
            end

        case 'AEDIM'
            varStruct.(varname).arcdependent=1;
            arcnum=retrieveodemodelinformation(ocStruct,'arcnum');
            arcidentifier=retrieveodemodelinformation(ocStruct,'arcidentifier');
            for ii=1:arcnum.value
                algebraicequationnum=retrieveodemodelinformation(ocStruct,'algebraicequationnum',arcidentifier.value{ii});
                varStruct.(varname).string{ii}=int2str(algebraicequationnum.value);
            end

        case 'DAEORDER'
            if ~isfield(varStruct,'AEDIM')
                varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'AEDIM');
            end
            varStruct.(varname).arcdependent=1;
            arcnum=retrieveodemodelinformation(ocStruct,'arcnum');
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            for ii=1:arcnum.value
                order=['[' int2str(ones(1,statenum.value))];
                aedim=str2double(varStruct.AEDIM.string{ii});
                if ~isnan(aedim) && aedim
                    order=[order '  ' int2str(zeros(1, aedim))];
                end
                order=[order ']'];
                varStruct.(varname).string{ii}=order;
            end

        case 'SUBSFLAG'
            varStruct.(varname).string='s';

        case 'EXOGENOUSFUNCTION'
            varStruct.(varname).string=getbasicname('exogenousfunction');

        case 'EXOGENOUSFUNCTIONTERM'
            try
                exogenousfunctionterm=algebraicodeterm2string(ocStruct,'exogenousfunctionterm',1,'');
                for jj=1:numel(exogenousfunctionterm.term)
                    if jj==1
                        varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                    else
                        varStruct.(varname).string{jj}=exogenousfunctionterm.term{jj};
                    end
                end
                varStruct.(varname).string{jj}=[varStruct.(varname).string{jj}];
                %end
                varStruct.(varname).multline=1;
                varStruct.(varname).vectorize=1;
                varStruct.(varname).arcdependent=0;
                varStruct.(varname).type='algterm';
            catch
                %lasterr
                varStruct.(varname).string=varname;
            end

        case 'SYMBOLICEXOGENOUSFUNCTIONTERM'
            symbolicexogenousfunctionterm=algebraicodeterm2string(ocStruct,'symbolicexogenousfunctionterm',0,'');
            for jj=1:numel(symbolicexogenousfunctionterm.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' symbolicexogenousfunctionterm.term{jj}];
                else
                    varStruct.(varname).string{jj}=symbolicexogenousfunctionterm.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';
            
            
        case 'SYMBOLICEXOGENOUSFUNCTIONNAMESUBSTITUTION'
            exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionname');
            exogenousfunctionstr=getbasicname('exogenousfunction');
            try
                varStruct.(varname).string=['out=subs(out,mystr2sym({''' exogenousfunctionname.value{1} '()''}),{' exogenousfunctionstr '});'];
            catch
                varStruct.(varname).string='';
            end
            
        case 'GENERATIONDATE'
            varStruct.(varname).string=datestr(now);

        case 'MODELNAME'
            varStruct.(varname).string=modelname(ocStruct);

        case 'MODELTYPE'
            varStruct.(varname).string=modeltype(ocStruct);

        case 'ARCIDENTIFIER'
            varStruct.(varname).arcdependent=1;
            arcidentifier=retrieveodemodelinformation(ocStruct,'arcidentifier');
            varStruct.(varname).string=arcidentifier.value;

        case 'NUMBEROFARCS'
            arcnum=retrieveodemodelinformation(ocStruct,'arcnum');
            varStruct.(varname).string=num2str(arcnum.value);

        case 'ARCARGUMENT'
            arcargument=retrieveodemodelinformation(ocStruct,'argument');
            varStruct.(varname).string=mat2str(arcargument.value);

        case 'EDGES'
            arcargument=retrieveodemodelinformation(ocStruct,'argument');
            arcnum=retrieveodemodelinformation(ocStruct,'arcnum');
            if arcnum.value>1
                nc=nchoosek(arcargument.value,2);
                varStruct.(varname).string=mat2str([nc;nc(:,[2 1])]);
            else
                varStruct.(varname).string='[]';
            end

        case 'ARCARGUMENT'
            arcargument=retrieveodemodelinformation(ocStruct,'argument');
            varStruct.(varname).string=mat2str(arcargument.value);

        case 'ISAUTONOMOUS'
            autonomousflag=retrieveodemodelinformation(ocStruct,'autonomous');
            varStruct.(varname).string=num2str(autonomousflag.value);

        case 'EXOGENOUSFUNCTIONNUM'
            exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
            varStruct.(varname).string=num2str(exogenousfunctionnum.value);

        case 'STATECOORD'
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=['1:' int2str(statenum.value)];

        case 'STATENUM'
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=int2str(statenum.value);

        case 'ZEROSONESTATENUM'
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=mat2str(zeros(1,statenum.value));

        case 'ZEROSTWOSTATENUM'
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            varStruct.(varname).string=mat2str(zeros(2,statenum.value));

        case 'DEFAULTPLOTCOMMAND'
            depvarname=getbasicname('dependent');
            varStruct.(varname).string=['h=plot(' depvarname '(1,:),' depvarname '(2,:));'];

        case 'DEFAULTPLOTTHREECOMMAND'
            depvarname=getbasicname('dependent');
            varStruct.(varname).string=['h=plot3(' depvarname '(1,:),' depvarname '(2,:),' depvarname '(3,:));'];


        case 'PARAMETERVALUES'
            % the strings used for assigning parametervalues from the par
            % argument in functions to the spcific paraemter variable names,
            % e.g. r=pararg(1); ...
            varStruct.(varname).multline=2;
            basename=getbasicname('parametervariables');
            parvar=algebraicodeterm2string(ocStruct,'parametervariables');
            for ii=1:numel(parvar.term)
                varStruct.(varname).string{ii}=[parvar.term{ii} '=' basename '(' num2str(ii) ');'];
            end

        case 'LATEXUSERDEPENDENTNAMEFIRST'
            statename=retrieveodemodelinformation(ocStruct,'statename');
            [val1 val2 val3]=regexp(statename.value{1},'([0-9]+)\>','split');
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' statename.value{1}(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end

        case 'LATEXUSERDEPENDENTNAMESECOND'
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            if statenum.value>1
                statename=retrieveodemodelinformation(ocStruct,'statename');
                labelname=statename.value{2};
            else
                labelname='';
            end
            [val1 val2 val3]=regexp(labelname,'([0-9]+)\>','split');
            if ~isempty(val2)
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '_' labelname(val2:val3) '$'];
            else
                varStruct.(varname).string=['$' latex(mystr2sym(val1{1})) '$'];
            end

        case 'LATEXUSERDEPENDENTNAMETHIRD'
            statenum=retrieveodemodelinformation(ocStruct,'statenum');
            if statenum.value<=2
                labelname='';
            else
                statename=retrieveodemodelinformation(ocStruct,'statename');
                labelname=statename.value{3};
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

        case 'SYMBOLICDYNAMICS'
            symbolicdynamics=algebraicodeterm2string(ocStruct,'symbolicdynamics',0,'');
            for jj=1:numel(symbolicdynamics.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' symbolicdynamics.term{jj}];
                else
                    varStruct.(varname).string{jj}=symbolicdynamics.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=0;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';

        case 'DYNAMICS'
            dynamics=algebraicodeterm2string(ocStruct,'dynamics',1,'');
            for jj=1:numel(dynamics.term)
                if jj==1
                    varStruct.(varname).string{jj}=['out=' dynamics.term{jj}];
                else
                    varStruct.(varname).string{jj}=dynamics.term{jj};
                end
            end
            varStruct.(varname).string{jj}=[varStruct.(varname).string{jj} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).vectorize=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';


        case 'DYNAMICSJACOBIAN'
            dynamicsjacobian=algebraicodeterm2string(ocStruct,'dynamicsjacobian');
            for ii=1:numel(dynamicsjacobian.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' dynamicsjacobian.term{ii}];
                else
                    varStruct.(varname).string{ii}=dynamicsjacobian.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).type='algterm';

        case 'DYNAMICSPARAMETERJACOBIAN'
            dynamicsparameterjacobian=algebraicodeterm2string(ocStruct,'dynamicsparameterjacobian');
            for ii=1:numel(dynamicsparameterjacobian.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' dynamicsparameterjacobian.term{ii}];
                else
                    varStruct.(varname).string{ii}=dynamicsparameterjacobian.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).type='algterm';


        case 'SYMBOLICDYNAMICSJACOBIAN'
            dynamicsjacobian=algebraicodeterm2string(ocStruct,'symbolicdynamicsjacobian');
            for ii=1:numel(dynamicsjacobian.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' dynamicsjacobian.term{ii}];
                else
                    varStruct.(varname).string{ii}=dynamicsjacobian.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';

        case 'SYMBOLICDYNAMICSPARAMETERJACOBIAN'
            symbolicdynamicsparameterjacobian=algebraicodeterm2string(ocStruct,'symbolicdynamicsparameterjacobian');
            for ii=1:numel(symbolicdynamicsparameterjacobian.term)
                if ii==1
                    varStruct.(varname).string{ii}=['out=' symbolicdynamicsparameterjacobian.term{ii}];
                else
                    varStruct.(varname).string{ii}=symbolicdynamicsparameterjacobian.term{ii};
                end
            end
            varStruct.(varname).string{ii}=[varStruct.(varname).string{ii} ';'];
            varStruct.(varname).multline=1;
            varStruct.(varname).arcdependent=0;
            varStruct.(varname).type='algterm';

        case 'DATAPATHNAME'
            varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct))) ''''];

        case 'RESULTFILE'
            varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct)),'SaveIntermediateResults') ''''];

        case 'GLOBALVARFILE'
            varStruct.(varname).string=['''' fullocmatfile(getocmatfolder('userdata',modeltype(ocStruct),modelname(ocStruct)),'SaveIntermediateResultsGlobalVariable') ''''];

        case 'INFODETAILS'
            varStruct.(varname).multline=2;
            varStruct.(varname).string={'%', ...
                ['% this file was automatically created: ' datestr(now)], ...
                '% written by Dieter Grass, 2013'};
    end
catch
    ocmatmsg('Problems with variable ''%s''.',varname)
    varStruct.(varname).string=varname;
end

