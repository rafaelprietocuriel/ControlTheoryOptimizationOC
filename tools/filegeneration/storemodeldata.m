function datafilename=storemodeldata(ocStruct,modeldatafolder)
%
%
fullfolder=fullocmatfile(modeldatafolder);
searchpath=path;
if isempty(strfind(searchpath,[fullfolder pathsep]))
    answer=input([strrep(fullfolder,filesep,'\\') ' is not on MATLAB path. Add it?  (y)/n: '],'s');
    if isempty(answer)
        % default value 'y'
        answer='y';
    end
    if strcmpi(answer,'n')
        disp(['Not added to MATLAB path.'])
    else
        addpath(fullfolder)
        savepath
    end
end
datafilename=fullfile(fullfolder,[modelname(ocStruct) 'ModelDataStructure.mat']);
save(datafilename,'ocStruct','-mat');