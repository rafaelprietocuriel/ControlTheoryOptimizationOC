function ocgTrjC=ocgradmultipath2cell(ocgMTrj)
ocgTrjC=[];
if isempty(ocgMTrj)
    return
end
var=variable(ocgMTrj);
var_fn=fieldnames(var);

switch octype(ocgMTrj)
    case 'static'
        d=multiplicity(ocgMTrj);
        ocgTrjC=cell(d,1);
        for ii=1:d
            ocgTrjStr=struct(ocgradtrajectory());
            ocgTrjStr.initargument=[];
            ocgTrjStr.argument=[];
            for jj=1:length(var_fn)
                ocgTrjStr.variable.(var_fn{jj})=var.(var_fn{jj})(:,ii);
            end
            ocgTrjStr.timehorizon=[];
            ocgTrjStr.modelname=modelname(ocgMTrj);
            ocgTrjStr.modelparameter=modelparameter(ocgMTrj);
            ocgTrjStr.solver=solver(ocgMTrj);
            ocgTrjC{ii}=ocgradtrajectory(ocgTrjStr);
        end
    otherwise
        multiplic=multiplicity(ocgMTrj);
        d=multiplic.number;
        ocgTrjC=cell(d,1);

        timestructure=multiplic.timestructure;

        initarg=initargument(ocgMTrj);
        initarg_fn=fieldnames(initarg);

        arg=argument(ocgMTrj);
        arg_fn=fieldnames(arg);
        T=timehorizon(ocgMTrj);

        for ii=1:d
            ocgTrjStr=struct(ocgradtrajectory());
            for jj=1:length(initarg_fn)
                ocgTrjStr.initargument.(initarg_fn{jj})=initarg.(initarg_fn{jj})(ii);
            end
            for jj=1:length(arg_fn)
                ocgTrjStr.argument.(arg_fn{jj})=arg.(arg_fn{jj})(:,timestructure(1,ii):timestructure(2,ii));
            end
            for jj=1:length(var_fn)
                ocgTrjStr.variable.(var_fn{jj})=var.(var_fn{jj})(:,timestructure(1,ii):timestructure(2,ii));
            end
            ocgTrjStr.timehorizon=T(ii);
            ocgTrjStr.modelname=modelname(ocgMTrj);
            ocgTrjStr.modelparameter=modelparameter(ocgMTrj);
            ocgTrjStr.solver=solver(ocgMTrj);
            ocgTrjC{ii}=ocgradtrajectory(ocgTrjStr);
        end
end