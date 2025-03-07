function odeObj=delete(odeObj,varargin)
%
% DELETE deletes specific results stored in ODEOBJ.
%
% ODEOBJ=DELETE(ODEOBJ,FIELDNAME,INDEX) deletes the results stored under FIELDNAME in the 
% ocResults, where INDEX is an array specifying the index of the results to be deleted.
%
% ODEOBJ=DELETE(ODEOBJ,FIELDNAME1,INDEX1,FIELDNAME2,INDEX2, ...) 

if isempty(odeObj)
    return
end
resultStruct=result(odeObj);
fieldname=fieldnames(resultStruct);
for ii=1:nargin-1
    if iscell(varargin{ii})
        deletecoord=varargin{ii+1};
        for jj=1:length(varargin{ii})
            deletefieldcoor(varargin{ii}{jj})
        end
    elseif ischar(varargin{ii})
        deletecoord=varargin{ii+1};
        deletefieldcoor(varargin{ii})
    end
end
    function deletefieldcoor(df)
        if ~strcmp(df,fieldname)
            warning([df ' is not an result.'])
            return
        else
            lr=length(resultStruct.(df));
            if ~isempty(deletecoord) && (any(max(deletecoord))>lr || any(min(deletecoord))<1 || any(rem(deletecoord,1)))
                error('Indices must be a positive integer vector.')
            end
            if isempty(deletecoord)
                deletecoord=1:lr;
            end
            resultStruct.(df)(deletecoord)=[];
        end
    end
odeObj.Result=resultStruct;


if ~nargout
    assignin('caller',inputname(1),odeObj);
else
    varargout{1}=odeObj;
end

end
