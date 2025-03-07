function transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,varname,overwrite,varargin)
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
    case 'statet'
        independentvar=getbasicname('independent');
        basename=getbasicname('dependent');
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii) '_' independentvar];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ',1)'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',1:end-1)'];
        end
    case 'statetp1'
        independentvar=getbasicname('independent');
        basename=getbasicname('dependent');
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii) '_' independentvar 'p1'];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ',2)'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',2:end)'];
        end
    case 'costatet'
        independentvar=getbasicname('independent');
        basename=getbasicname('dependent');
        costatenum=retrievediffmodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii+ocStruct.variable.state.num;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord) '_' independentvar];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ',1)'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',1:end-1)'];
        end
    case 'costatetp1'
        independentvar=getbasicname('independent');
        basename=getbasicname('dependent');
        costatenum=retrievediffmodelinformation(ocStruct,'costatenum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        for ii=1:costatenum.value
            costatecoord=ii+ocStruct.variable.state.num;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord) '_' independentvar 'p1'];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(costatecoord) ',2)'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(costatecoord) ',2:end)'];
        end
    case 'controltp1'
        independentvar=getbasicname('independent');
        basename=getbasicname('control');
        implicitbasename=getbasicname('dependent');
        controlnum=retrievediffmodelinformation(ocStruct,'controlnum');
        arcnum=retrievediffmodelinformation(ocStruct,'arcnum');
        statenum=retrievediffmodelinformation(ocStruct,'statenum');
        arcid=retrievediffmodelinformation(ocStruct,'arcidentifier');
        transStruct.(varname).arcdependent=1;
        transStruct.(varname).basicname=basename;
        for ii=1:controlnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii) '_' independentvar 'p1'];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ',1)'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcid.value{ii});
            implicitnonlinearcontrolindex=retrievediffmodelinformation(ocStruct,'implicitnonlinearcontrolindex',arcid.value{ii});
            implicitcounter=0;
            for jj=1:controlnum.value
                if any(implicitnonlinearcontrolindex.value==jj)
                    implicitcounter=implicitcounter+1;
                    controlcoord=2*statenum.value+implicitcounter;
                    transStruct.(varname).(arc).name{jj}=[implicitbasename int2str(controlcoord)];
                    transStruct.(varname).(arc).arrayname{jj}=[implicitbasename '(' int2str(controlcoord) ')'];
                    transStruct.(varname).(arc).vectorizedname{jj}=[implicitbasename '(' int2str(controlcoord) ',:)'];
                else
                    transStruct.(varname).(arc).name{jj}=[basename int2str(jj) '_' independentvar 'p1'];
                    transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(jj) ',1)'];
                    transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(jj) ',1:end-1)'];
                end
            end
        end
    case 'lagrangemultcctp1'
        independentvar=getbasicname('independent');
        cconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if ~cconstraintnum.value
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
            return
        end
        basename=getbasicname('lagrangemultcc');
        transStruct.(varname).basicname=basename;
        for ii=1:cconstraintnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii) '_' independentvar 'p1'];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ',1)'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
    case 'exogenousfunction'
        basename=getbasicname('exogenousfunction');
        exogenousfunctionnum=retrievediffmodelinformation(ocStruct,'exogenousfunctionnum');
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
    case 'vectorizedsize'
        transStruct.(varname).string=['size(' getbasicname('dependent') ',2)-1'];
        transStruct.(varname).arcdependent=0;
    case 'composedvar'
        transStruct.(varname).arcdependent=1;
        if ~isfield(transStruct,'statet')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'statet');
        end
        if ~isfield(transStruct,'costatet')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'costatet');
        end
        if ~isfield(transStruct,'statetp1')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'statetp1');
        end
        if ~isfield(transStruct,'costatetp1')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'costatetp1');
        end
        if ~isfield(transStruct,'controltp1')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'controltp1');
        end
        if ~isfield(transStruct,'lagrangemultcctp1')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'lagrangemultcctp1');
        end
        if ~isfield(transStruct,'exogenousfunction')
            transStruct=generatetransformationstruct4standarddiffmodel(ocStruct,transStruct,'exogenousfunction');
        end
        for ii=1:ocStruct.arc.num
            if ii==1
                parametername=retrievediffmodelinformation(ocStruct,'parametername');
                parameternum=retrievediffmodelinformation(ocStruct,'parameternum');
                cconstraintnum=retrievediffmodelinformation(ocStruct,'inequalitycontrolconstraintnum');
                statenamet=retrievediffmodelinformation(ocStruct,'statenamet');
                costatenamet=retrievediffmodelinformation(ocStruct,'costatenamet');
                statenametp1=retrievediffmodelinformation(ocStruct,'statenametp1');
                costatenametp1=retrievediffmodelinformation(ocStruct,'costatenametp1');
                controlnametp1=retrievediffmodelinformation(ocStruct,'controlnametp1');
                lagrangemultccnametp1=retrievediffmodelinformation(ocStruct,'lagrangemultipliercontrolnametp1');
                exogenousfunctionnum=retrievediffmodelinformation(ocStruct,'exogenousfunctionnum');
                exogenousfunctionname=retrievediffmodelinformation(ocStruct,'exogenousfunctionname');
                transStruct.(varname).parametername=parametername.value;
                for jj=1:parameternum.value
                    transStruct.(varname).vectorizedparametername{jj}=['repmat(' parametername.value{jj} ',1,size(' transStruct.statet.basicname ',2)-1)'];
                end
                transStruct.(varname).username=[statenamet.value costatenamet.value statenametp1.value costatenametp1.value controlnametp1.value];
                transStruct.(varname).standardname=[transStruct.statet.name transStruct.costatet.name transStruct.statetp1.name transStruct.costatetp1.name transStruct.controltp1.name];
                transStruct.(varname).standardarrayname=[transStruct.statet.arrayname transStruct.costatet.arrayname transStruct.statetp1.arrayname transStruct.costatetp1.arrayname transStruct.controltp1.arrayname];
                transStruct.(varname).standardvectorizedname=[transStruct.statet.vectorizedname transStruct.costatet.vectorizedname transStruct.statetp1.vectorizedname transStruct.costatetp1.vectorizedname transStruct.controltp1.vectorizedname];
                if cconstraintnum.value
                    transStruct.(varname).username=[transStruct.(varname).username lagrangemultccnametp1.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.lagrangemultcctp1.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.lagrangemultcctp1.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.lagrangemultcctp1.vectorizedname];
                end
                if exogenousfunctionnum.value
                    transStruct.(varname).username=[transStruct.(varname).username exogenousfunctionname.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.exogenousfunction.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.exogenousfunction.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.exogenousfunction.vectorizedname];
                end
            end
            arc=arcidentifier2field(ocStruct.arc.identifier{ii});
            transStruct.(varname).(arc).name=[transStruct.statet.name transStruct.costatet.name transStruct.statetp1.name transStruct.costatetp1.name transStruct.controltp1.(arc).name transStruct.lagrangemultcctp1.name transStruct.exogenousfunction.name];
            transStruct.(varname).(arc).arrayname=[transStruct.statet.arrayname transStruct.costatet.arrayname transStruct.statetp1.arrayname transStruct.costatetp1.arrayname transStruct.controltp1.(arc).arrayname transStruct.lagrangemultcctp1.arrayname transStruct.exogenousfunction.arrayname];
            transStruct.(varname).(arc).vectorizedname=[transStruct.statet.vectorizedname transStruct.costatet.vectorizedname transStruct.statetp1.vectorizedname transStruct.costatetp1.vectorizedname transStruct.controltp1.(arc).vectorizedname transStruct.lagrangemultcctp1.vectorizedname transStruct.exogenousfunction.vectorizedname];
        end
end        