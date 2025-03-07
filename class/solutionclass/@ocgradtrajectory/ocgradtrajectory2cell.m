function ocgTrjC=ocgradtrajectory2cell(ocgMTrj)
ocgTrjC=[];
if isempty(ocgMTrj)
    return
end
d=multiplicity(ocgMTrj);
ocgTrjC=cell(d,1);

for ii=1:d
    var=variable(ocgMTrj(ii));
    var_fn=fieldnames(var);

    switch octype(ocgMTrj(ii))
        case 'static'
            ocgTrjStr=struct(ocgradtrajectory());
            ocgTrjStr.initargument=[];
            ocgTrjStr.argument=[];
            for jj=1:length(var_fn)
                ocgTrjStr.variable.(var_fn{jj})=var.(var_fn{jj})(:,ii);
            end
            ocgTrjStr.timehorizon=[];
            ocgTrjStr.modelname=modelname(ocgMTrj(ii));
            ocgTrjStr.modelparameter=modelparameter(ocgMTrj(ii));
            ocgTrjStr.solver=solver(ocgMTrj(ii));
            ocgTrjC{ii}=ocgradtrajectory(ocgTrjStr);
        otherwise
            ocgTrjC{ii}=ocgradtrajectory(ocgMTrj(ii));
    end
end