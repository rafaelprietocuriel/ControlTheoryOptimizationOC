function b = fieldcontrol(fieldname,varargin)

if nargin == 2
    value = varargin{1};
    flag=1;
else
    flag=0;
end

if ~isstr(fieldname)
    error('Fieldname must be string.')
end



switch fieldname
    case 'AutoDataFileGeneration'
        if flag
            if ~strcmp(value,'y') && ~strcmp(value,'n') && ~strncmp(value,'a')
                error('Incorrect option value for field ''AutoDataFileGeneration''.')
            end
        end
    case 'OutputFormatODE'
        if flag
            if ~strcmp(value,'ocmat') && ~strcmp(value,'matlab')
                error('Incorrect option value for field ''OutputFormatODE''.')
            end
        end
    otherwise
        error('Incorrect fieldname.')
end
b = 1;
        