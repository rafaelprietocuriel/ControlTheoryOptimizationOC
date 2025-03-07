function varargout=removeresult(ppdeObj,fieldname,varargin)
%
% REMOVERESULT removes all results from a Result field.
%
% OCOBJ=REMOVERESULT(OCOBJ,FIELDNAME) removes the field FIELDNAME from the
% Result.
%
% OCOBJ=REMOVERESULT(OCOBJ,FIELDNAME,'FORCE') forces to remove the field
% without a warning.

force=[];

if isempty(ppdeObj)
    return
end
resultStruct=result(ppdeObj);
if ~isfield(resultStruct,fieldname)
    warning([fieldname ' is not a field of the Result.'])
    return
end

if nargin==2
    force=0;
end

if nargin==3
    if ischar(varargin{1})
        force=strcmpi(varargin{1}(1),'f');
    else
        force=varargin{1};
    end
end

if ~force
    answer=input(['Do you really want to remove field ''' fieldname ''' from the Result?  y/(n): '],'s');
    if isempty(answer)
        % default value 'n'
        answer='n';
    end
    if strcmpi(answer,'n')
        disp('Return without removing.')
        return
    end
end
ppdeObj.Result=rmfield(ppdeObj.Result,fieldname);

if ~nargout
    assignin('caller',inputname(1),ppdeObj);
else
    varargout{1}=ppdeObj;
end

end
