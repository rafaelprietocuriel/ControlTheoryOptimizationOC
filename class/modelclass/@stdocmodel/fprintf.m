function varargout=fprintf(ocObj,filename,varargin)

modelpath=fullocmatfile(getocmatfolder('specificmodel','',modelname(ocObj)));
if isempty(filename)
   return 
end

absolutefilename=fullfile(modelpath,[modelname(ocObj) filename '.m']);

if ~exist(absolutefilename,'file')
    [fidwrite,msg]=fopen(absolutefilename,'w+');
    if fidwrite==-1
        ocmatmsg(msg)
        return
    end
    varStruct=generatevariablestruct4standardmodel(ocObj.Model,[],'PARAMETERVALUES',1);
    fprintf(fidwrite, ['function out=' modelname(ocObj) filename '( )\n']);
    fprintf(fidwrite, '%s\t%s:\n\n',repmat('%',1,2),datestr(date,'dd-mmm-yyyy'));
    for ii=1:numel(varStruct.PARAMETERVALUES.string)
        fprintf(fidwrite, '%s\n',varStruct.PARAMETERVALUES.string{ii});
    end
    fprintf(fidwrite, '\n');
else
    [fidwrite,msg]=fopen(absolutefilename,'a+');
    if fidwrite==-1
        ocmatmsg(msg)
        return
    end
end

count=fprintf(fidwrite,varargin{:});
fclose(fidwrite);

if nargout==1
    varargout{1}=count;
end