function newfilename=makefilenamefromtemplate(templatename,modelname)

newfilename=fullocmatfile(fullfile('standardmodel','out',[strrep(templatename,'standardmodel',modelname) '.m']));
