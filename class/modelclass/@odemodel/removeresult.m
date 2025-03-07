function varargout=removeresult(odeObj,fieldname,varargin)
%
% REMOVERESULT removes all results from a Result field.
%
% ODEOBJ=REMOVERESULT(ODEOBJ,FIELDNAME) removes the field FIELDNAME from the
% Result.
%
% ODEOBJ=REMOVERESULT(ODEOBJ,FIELDNAME,'FORCE') forces to remove the field
% without a warning.

force=[];

if isempty(odeObj)
    return
end
resultStruct=result(odeObj);
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
odeObj.Result=rmfield(odeObj.Result,fieldname);

if ~nargout
    assignin('caller',inputname(1),odeObj);
else
    varargout{1}=odeObj;
end

end