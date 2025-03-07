function transStruct=generatetransformationstruct4odemodel(ocStruct,transStruct,varname,overwrite,varargin)
%
% in this structure the information for the transformation of optimal
% control specific variables to default names, which are used in the
% automatically generated files are stored.
% 
if nargin<=3
    overwrite=0;
end

if ~overwrite
    if isfield(transStruct,varname)
        return
    end
end

switch varname
    case 'independent'
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).name=getbasicname(varname);
    case 'state'
        basename=getbasicname('dependent');
        statenum=retrieveodemodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end

    case 'composedvar'
        transStruct.(varname).arcdependent=1;
        if ~isfield(transStruct,'state')
            transStruct=generatetransformationstruct4odemodel(ocStruct,transStruct,'state');
        end
        if ~isfield(transStruct,'exogenousfunction')
            transStruct=generatetransformationstruct4odemodel(ocStruct,transStruct,'exogenousfunction');
        end
        for ii=1:ocStruct.arc.num
            if ii==1
                parametername=retrieveodemodelinformation(ocStruct,'parametername');
                parameternum=retrieveodemodelinformation(ocStruct,'parameternum');
                statename=retrieveodemodelinformation(ocStruct,'statename');
                exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
                exogenousfunctionname=retrieveodemodelinformation(ocStruct,'exogenousfunctionnamewithargument');
                transStruct.(varname).parametername=parametername.value;
                for jj=1:parameternum.value
                    transStruct.(varname).vectorizedparametername{jj}=['repmat(' parametername.value{jj} ',1,size(' transStruct.state.basicname ',2))'];
                end
                transStruct.(varname).username=[statename.value];
                transStruct.(varname).standardname=[transStruct.state.name];
                transStruct.(varname).standardarrayname=[transStruct.state.arrayname];
                transStruct.(varname).standardvectorizedname=[transStruct.state.vectorizedname];
                if exogenousfunctionnum.value
                    transStruct.(varname).username=[transStruct.(varname).username exogenousfunctionname.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.exogenousfunction.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.exogenousfunction.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.exogenousfunction.vectorizedname];
                end
            end
            arc=arcidentifier2field(ocStruct.arc.identifier{ii});
            transStruct.(varname).(arc).name=[transStruct.state.name transStruct.exogenousfunction.name];
            transStruct.(varname).(arc).arrayname=[transStruct.state.arrayname transStruct.exogenousfunction.arrayname];
            transStruct.(varname).(arc).vectorizedname=[transStruct.state.vectorizedname transStruct.exogenousfunction.vectorizedname];
        end
        
    case 'exogenousfunction'
        basename=getbasicname('exogenousfunction');
        exogenousfunctionnum=retrieveodemodelinformation(ocStruct,'exogenousfunctionnum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        if exogenousfunctionnum.value
            for ii=1:exogenousfunctionnum.value
                transStruct.(varname).name{ii}=[basename int2str(ii)];
                transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
            end
        else
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
        end

end        