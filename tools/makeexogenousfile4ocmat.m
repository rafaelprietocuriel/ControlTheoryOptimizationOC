function modelfiles=makeexogenousfile4ocmat(modelname,odes,moveflag)
%
% MAKEEXOGENOUSFILE4OCMAT generates files for exogenous dynamics:
% modelnameExogenousDynamics: containing the dynamics
% modelnameExogenousJacobian: containing the Jacobian of the dynamics with
% respect to the states, costates and exogenous state names
% modelnameExogenousParameterJacobian: containing the Jacobian of the
% dynamics with respect to the parameter names
%
% MAKEEXOGENOUSFILE4OCMAT(OCOBJ,ODES) standardmodel OCOBJ, symobolid ODES in the
% form sym('Dx=f(x)')
%
% MAKEEXOGENOUSFILE4OCMAT(OCOBJ,ODES,MOVEFLAG) MOVEFLAG=1 the files are moved to
% the model folder, MOVEFLAG=0 are stored in standard output folder. Default
% value is MOVEFLAG=1.
%
% Example usage:
% odes=sym('[DN=N*(b*l-delta),Dc=x^(sigma/beta)*h^(sigma*mu/beta+1)*(1-l-e)*N]')
% m=stdocmodel('skibajonesfinal');
% modelfiles=makeexogenousfile4ocmat(m,odes)

if nargin==2
    moveflag=1;
end
ocStruct=loadmodeldata(modelname);
ocObj=stdocmodel(modelname,[],[],0);

varStruct=generatevariablestruct4standardmodel(ocStruct,[],'INFODETAILS');
varStruct=generatevariablestruct4standardmodel(ocStruct,varStruct,'PARAMETERVALUES');

outdir=fullocmatfile(getocmatfolder('out',modeltype(ocObj)));

modelfiles(1).name='ExogenousDynamics';
modelfiles(2).name='ExogenousJacobian';
modelfiles(3).name='ExogenousParameterJacobian';
modelfiles(4).name='ExogenousDynamicsDerivativeTime';

arcn=arcnum(ocObj);

x=state(ocObj);
l=costate(ocObj);

parn=parametername(ocObj);

numodes=length(odes);

exogenousname=cell(1,numodes);
dynamics=cell(1,numodes);
jacob=cell(1,arcn);
parjacob=cell(1,arcn);

for ii=1:numodes
    [ls,dynamics{ii},exogenousname{ii}]=splitode(char(odes(ii)));
end
X=[x(:);l(:);exogenousname(:)];
ocStruct.variable.exogenousstate.name=exogenousname;
ocStruct.variable.exogenousstate.num=numodes;

optdynamics=cell(1,arcn);
optdynamicsstr=cell(arcn,numodes);
jacobstr=cell(arcn,numodes);
parjacobstr=cell(arcn,numodes);
timeder=cell(1,arcn);
timederstr=cell(arcn,numodes);

for ii=1:arcn
    arcid=ii-1;
    optdynamics{ii}=term2optterm_standardmodel(ocStruct,dynamics,arcid);
    jacob{ii}=sym(removematrixstring(ocmatjacobian(cell2vectorstring(optdynamics{ii}),cell2vectorstring(X))));
    parjacob{ii}=sym(removematrixstring(ocmatjacobian(cell2vectorstring(optdynamics{ii}),cell2vectorstring(parn))));
    timeder{ii}=sym(removematrixstring(ocmatjacobian(cell2vectorstring(optdynamics{ii}),'[t]')));
    for jj=1:numodes
        transterm=term2standardvariable(ocStruct,jacob{ii}(jj,:),arcid,2,1);
        jacobstr{ii,jj}=transterm{1}(2:end-1);
        transterm=term2standardvariable(ocStruct,parjacob{ii}(jj,:),arcid,2,1);
        parjacobstr{ii,jj}=transterm{1}(2:end-1);
        transterm=term2standardvariable(ocStruct,optdynamics{ii}{jj},arcid,2,1);
        optdynamicsstr{ii,jj}=transterm{1};
        transterm=term2standardvariable(ocStruct,timeder{ii}(jj),arcid,2,1);
        timederstr{ii,jj}=transterm{1};
    end
end
modelfiles(1).string=optdynamicsstr;
modelfiles(2).string=jacobstr;
modelfiles(3).string=parjacobstr;
modelfiles(4).string=timederstr;

for ii=1:length(modelfiles)
    modelfiles(ii).filename=makefile(modelfiles(ii).name,modelfiles(ii).string);
end

if moveflag
    moveocmatfiles(ocStruct,modelfiles);
end

    function totfilename=makefile(funcname,funcstr)
        filename=[modelname funcname];
        totfilename=fullfile(outdir,[filename '.m']);
        fid=fopen(totfilename,'w');
        fn=fieldnames(varStruct);
        
        % Header
        fprintf(fid,'function out=%s(t,depvar,par,arcid)\n',filename);
        for kk=1:length(fn)
            if ~varStruct.(fn{kk}).arcdependent
                for ll=1:length(varStruct.(fn{kk}).string)
                    fprintf(fid,'%s\n',varStruct.(fn{kk}).string{ll});
                end
            end
        end
        fprintf(fid,'\nswitch arcid\n');
        for kk=0:arcn-1
            fprintf(fid,'\tcase %s\n',num2str(kk));
            fprintf(fid,'\t\tout=');
            for ll=1:numodes
                if isempty(strfind(funcstr{kk+1,ll},'depvar')) && isempty(regexp(funcstr{kk+1,ll},'\<t\>','ONCE'))
                    funcstr{kk+1,ll}=sprintf('repmat(%s,1,length(t))',funcstr{kk+1,ll});
                end
                if ll==1
                    if numodes>1
                        fprintf(fid,'[%s; ...\n',funcstr{kk+1,ll});
                    else
                        fprintf(fid,'%s;\n',funcstr{kk+1,ll});
                    end
                elseif ll<numodes
                    fprintf(fid,'\t\t\t%s; ...\n',funcstr{kk+1,ll});
                else
                    fprintf(fid,'\t\t\t%s];\n',funcstr{kk+1,ll});
                end
            end
        end
        fprintf(fid,'end');
        fclose(fid);
    end
    function [ls,rs,varname]=splitode(ode)
        split=regexp(ode,'=','split');
        ls=split{1};
        rs=split{2};

        % remove D
        varname=ls(2:end);
    end

end