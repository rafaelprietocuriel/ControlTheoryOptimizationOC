function transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,varname,overwrite,varargin)
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
    case 'space'
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).name=getbasicname(varname);
    case 'time'
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).name=getbasicname(varname);
    case 'independent'
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).name=getbasicname(varname);
    case 'state'
        basename=getbasicname('dependent');
        gridnum=getbasicname('femdatagridnum');
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        statedependence=retrieveppdemodelinformation(ocStruct,'statedependence');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        ctr=0;
        for ii=1:statenum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            if strcmp(statedependence.value{ii},'2')
                if ii==1
                    transStruct.(varname).arrayname{ii}=[basename '(1:' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(1:' gridnum ',:)'];
                elseif ii==2
                    transStruct.(varname).arrayname{ii}=[basename '(' gridnum '+1:' int2str(ii) '*' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' gridnum '+1:' int2str(ii) '*' gridnum ',:)'];
                else
                    transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum ',:)'];
                end
            elseif strcmp(statedependence.value{ii},'0')
                ctr=ctr+1;
                base=[int2str(ii-1) '*' gridnum];
                transStruct.(varname).arrayname{ii}=[basename '(' base '+' int2str(ctr) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' base '+' int2str(ctr) ',:)'];
            end
        end
    case 'costate'
        basename=getbasicname('dependent');
        costatenum=retrieveppdemodelinformation(ocStruct,'costatenum');
        statenum=retrieveppdemodelinformation(ocStruct,'statenum');
        statedependence=retrieveppdemodelinformation(ocStruct,'statedependence');
        baseidx0=sum(strcmp(statedependence.value,'0'));
        baseidx2=sum(strcmp(statedependence.value,'2'));
        costatedependence=retrieveppdemodelinformation(ocStruct,'costatedependence');
        gridnum=getbasicname('femdatagridnum');
        transStruct.(varname).arcdependent=0;
        transStruct.(varname).basicname=basename;
        ctr=0;
        for ii=1:costatenum.value
            costatecoord=ii+statenum.value;
            transStruct.(varname).name{ii}=[basename int2str(costatecoord)];
            if strcmp(costatedependence.value{ii},'2')
                if baseidx0>0
                    base=[int2str(baseidx2) '*' gridnum '+' int2str(baseidx0)];
                else
                    base=[int2str(baseidx2) '*' gridnum];
                end
                if ii==1
                    transStruct.(varname).arrayname{ii}=[basename '(' base '+(1:' gridnum '))'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' base '+(1:' gridnum '),:)'];
                elseif ii==2
                    transStruct.(varname).arrayname{ii}=[basename '(' base '+(' gridnum '+1:' int2str(ii) '*' gridnum '))'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' base '+(' gridnum '+1:' int2str(ii) '*' gridnum '),:)'];
                else
                    transStruct.(varname).arrayname{ii}=[basename '(' base '+(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum '))'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' base '+(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum '),:)'];
                end
            elseif strcmp(costatedependence.value{ii},'0')
                ctr=ctr+1;
                base=[int2str(ii-1) '*' gridnum];
                transStruct.(varname).arrayname{ii}=[basename '(' base '+' int2str(ctr) ')'];
                transStruct.(varname).vectorizedname{ii}=[basename '(' base '+' int2str(ctr) ',:)'];
            end
        end
    case 'control'
        basename=getbasicname(varname);
        controlnum=retrieveppdemodelinformation(ocStruct,'controlnum');
        controldependence=retrieveppdemodelinformation(ocStruct,'controldependence');
        arcnum=retrieveppdemodelinformation(ocStruct,'arcnum');
        arcid=retrieveppdemodelinformation(ocStruct,'arcidentifier');
        gridnum=getbasicname('femdatagridnum');
        transStruct.(varname).arcdependent=1;
        transStruct.(varname).basicname=basename;
        ctr=0;
        for ii=1:controlnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            if strcmp(controldependence.value{ii},'2')
                if ii==1
                    transStruct.(varname).arrayname{ii}=[basename '(1:' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(1:' gridnum ',:)'];
                elseif ii==2
                    transStruct.(varname).arrayname{ii}=[basename '(' gridnum '+1:' int2str(ii) '*' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' gridnum '+1:' int2str(ii) '*' gridnum ',:)'];
                else
                    transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii-1) '*' gridnum '+1:' int2str(ii) '*' gridnum ',:)'];
                end
            elseif strcmp(controldependence.value{ii},'0')
                ctr=ctr+1;
                if ii==1
                    transStruct.(varname).arrayname{ii}=[basename '(1:' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(1:' gridnum ',:)'];
                elseif ii==2
                    transStruct.(varname).arrayname{ii}=[basename '(' gridnum '+1:' int2str(ctr) '*' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' gridnum '+1:' int2str(ctr) '*' gridnum ',:)'];
                else
                    transStruct.(varname).arrayname{ii}=[basename '(' int2str(ctr-1) '*' gridnum '+1:' int2str(ctr) '*' gridnum ')'];
                    transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ctr-1) '*' gridnum '+1:' int2str(ctr) '*' gridnum ',:)'];
                end
                %                 if ii==1
                %                     transStruct.(varname).arrayname{ii}=[basename '(' int2str(ctr) ')'];
                %                     transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ctr) ',:)'];
                %                 elseif ii==2
                %                     transStruct.(varname).arrayname{ii}=[basename '(' gridnum '+' int2str(ctr) ')'];
                %                     transStruct.(varname).vectorizedname{ii}=[basename '(' gridnum '+' int2str(ctr) ',:)'];
                %                 else
                %                     base=[int2str(ii-1) '*' gridnum];
                %                     transStruct.(varname).arrayname{ii}=[basename '(' base '+' int2str(ctr) ')'];
                %                     transStruct.(varname).vectorizedname{ii}=[basename '(' base '+' int2str(ctr) ',:)'];
                %                 end
            end
        end
        ctr=0;
        for ii=1:arcnum.value
            arc=arcidentifier2field(arcid.value{ii});
            for jj=1:controlnum.value
                transStruct.(varname).(arc).name{jj}=[basename int2str(jj)];
                if strcmp(controldependence.value{jj},'2')
                    if jj==1
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(1:' gridnum ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(1:' gridnum ',:)'];
                    elseif jj==2
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(' gridnum '+1:' int2str(jj) '*' gridnum ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' gridnum '+1:' int2str(jj) '*' gridnum ',:)'];
                    else
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(jj-1) '*' gridnum '+1:' int2str(jj) '*' gridnum ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(jj-1) '*' gridnum '+1:' int2str(jj) '*' gridnum ',:)'];
                    end
                elseif strcmp(controldependence.value{jj},'0')
                    ctr=ctr+1;
                    if jj==1
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(1:' gridnum ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(1:' gridnum ',:)'];
                    elseif jj==2
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(' gridnum '+1:' int2str(ctr) '*' gridnum ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' gridnum '+1:' int2str(ctr) '*' gridnum ',:)'];
                    else
                        transStruct.(varname).(arc).arrayname{jj}=[basename '(' int2str(ctr-1) '*' gridnum '+1:' int2str(ctr) '*' gridnum ')'];
                        transStruct.(varname).(arc).vectorizedname{jj}=[basename '(' int2str(ctr-1) '*' gridnum '+1:' int2str(ctr) '*' gridnum ',:)'];
                    end
                end
            end
        end
        ctr=0;
        for jj=1:controlnum.value
            transStruct.(varname).stdarc.name{jj}=[basename int2str(jj)];
            if strcmp(controldependence.value{jj},'2')
                if jj==1
                    transStruct.(varname).stdarc.arrayname{jj}=[basename '(1:' gridnum ')'];
                    transStruct.(varname).stdarc.vectorizedname{jj}=[basename '(1:' gridnum ',:)'];
                elseif jj==2
                    transStruct.(varname).stdarc.arrayname{jj}=[basename '(' gridnum '+1:' int2str(jj) '*' gridnum ')'];
                    transStruct.(varname).stdarc.vectorizedname{jj}=[basename '(' gridnum '+1:' int2str(jj) '*' gridnum ',:)'];
                else
                    transStruct.(varname).stdarc.arrayname{jj}=[basename '(' int2str(jj-1) '*' gridnum '+1:' int2str(jj) '*' gridnum ')'];
                    transStruct.(varname).stdarc.vectorizedname{jj}=[basename '(' int2str(jj-1) '*' gridnum '+1:' int2str(jj) '*' gridnum ',:)'];
                end
            elseif strcmp(controldependence.value{jj},'0')
                    ctr=ctr+1;
                    if jj==1
                        transStruct.(varname).stdarc.arrayname{jj}=[basename '(1:' gridnum ')'];
                        transStruct.(varname).stdarc.vectorizedname{ii}=[basename '(1:' gridnum ',:)'];
                    elseif jj==2
                        transStruct.(varname).stdarc.arrayname{jj}=[basename '(' gridnum '+1:' int2str(ctr) '*' gridnum ')'];
                        transStruct.(varname).stdarc.vectorizedname{jj}=[basename '(' gridnum '+1:' int2str(ctr) '*' gridnum ',:)'];
                    else
                        transStruct.(varname).stdarc.arrayname{jj}=[basename '(' int2str(ctr-1) '*' gridnum '+1:' int2str(ctr) '*' gridnum ')'];
                        transStruct.(varname).stdarc.vectorizedname{jj}=[basename '(' int2str(ctr-1) '*' gridnum '+1:' int2str(ctr) '*' gridnum ',:)'];
                    end
            end
        end
        
    case 'lagrangemultcc'
        cconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
        if ~cconstraintnum.value
            transStruct.(varname).name=[];
            transStruct.(varname).arrayname=[];
            transStruct.(varname).vectorizedname=[];
            return
        end
        basename=getbasicname(varname);
        transStruct.(varname).basicname=basename;
        for ii=1:cconstraintnum.value
            transStruct.(varname).name{ii}=[basename int2str(ii)];
            transStruct.(varname).arrayname{ii}=[basename '(' int2str(ii) ')'];
            transStruct.(varname).vectorizedname{ii}=[basename '(' int2str(ii) ',:)'];
        end
        
    case 'vectorizedsize'
        basename=getbasicname('dependent');
        transStruct.(varname).string=['size(' basename ',2)'];
        transStruct.(varname).arcdependent=0;
        
    case 'composedvar'
        transStruct.(varname).arcdependent=1;
        if ~isfield(transStruct,'state')
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'state');
        end
        if ~isfield(transStruct,'costate')
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'costate');
        end
        if ~isfield(transStruct,'control')
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'control');
        end
        if ~isfield(transStruct,'lagrangemultcc')
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'lagrangemultcc');
        end
        if ~isfield(transStruct,'lagrangemultsc')
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'lagrangemultsc');
        end
        if ~isfield(transStruct,'exogenousfunction')
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'exogenousfunction');
            transStruct=generatetransformationstruct4ppdemodel(ocStruct,transStruct,'exogenousfunctiondx');
        end
        for ii=1:ocStruct.arc.num
            if ii==1
                parametername=retrieveppdemodelinformation(ocStruct,'parametername');
                parameternum=retrieveppdemodelinformation(ocStruct,'parameternum');
                cconstraintnum=retrieveppdemodelinformation(ocStruct,'inequalitycontrolconstraintnum');
                statename=retrieveppdemodelinformation(ocStruct,'statename');
                costatename=retrieveppdemodelinformation(ocStruct,'costatename');
                controlname=retrieveppdemodelinformation(ocStruct,'controlname');
                lagrangemultccname=retrieveppdemodelinformation(ocStruct,'lagrangemultipliercontrolname');
                exogenousfunctionnum=retrieveppdemodelinformation(ocStruct,'exogenousfunctionnum');
                exogenousfunctionname=retrieveppdemodelinformation(ocStruct,'exogenousfunctionname');
                transStruct.(varname).parametername=parametername.value;
                for jj=1:parameternum.value
                    transStruct.(varname).vectorizedparametername{jj}=['repmat(' parametername.value{jj} ',1,size(' transStruct.state.basicname ',2))'];
                end
                transStruct.(varname).username=[statename.value costatename.value controlname.value];
                transStruct.(varname).standardname=[transStruct.state.name transStruct.costate.name transStruct.control.name];
                transStruct.(varname).standardarrayname=[transStruct.state.arrayname transStruct.costate.arrayname transStruct.control.arrayname];
                transStruct.(varname).standardvectorizedname=[transStruct.state.vectorizedname transStruct.costate.vectorizedname transStruct.control.vectorizedname];
                if cconstraintnum.value
                    transStruct.(varname).username=[transStruct.(varname).username lagrangemultccname.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.lagrangemultcc.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.lagrangemultcc.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.lagrangemultcc.vectorizedname];
                end
                if exogenousfunctionnum.value
                    transStruct.(varname).username=[transStruct.(varname).username exogenousfunctionname.value];
                    transStruct.(varname).standardname=[transStruct.(varname).standardname transStruct.exogenousfunction.name];
                    transStruct.(varname).standardarrayname=[transStruct.(varname).standardarrayname transStruct.exogenousfunction.arrayname];
                    transStruct.(varname).standardvectorizedname=[transStruct.(varname).standardvectorizedname transStruct.exogenousfunction.vectorizedname];
                end
            end
            arc=arcidentifier2field(ocStruct.arc.identifier{ii});
            transStruct.(varname).(arc).name=[transStruct.state.name transStruct.costate.name transStruct.control.(arc).name transStruct.lagrangemultcc.name transStruct.exogenousfunction.name];
            transStruct.(varname).(arc).arrayname=[transStruct.state.arrayname transStruct.costate.arrayname transStruct.control.(arc).arrayname  transStruct.lagrangemultcc.arrayname  transStruct.exogenousfunction.arrayname];
            transStruct.(varname).(arc).vectorizedname=[transStruct.state.vectorizedname transStruct.costate.vectorizedname transStruct.control.(arc).vectorizedname  transStruct.lagrangemultcc.vectorizedname  transStruct.exogenousfunction.vectorizedname];
        end
    case 'exogenousfunction'
        basename=getbasicname('exogenousfunction');
        exogenousfunctionnum=retrieveppdemodelinformation(ocStruct,'exogenousfunctionnum');
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