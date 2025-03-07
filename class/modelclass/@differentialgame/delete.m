function varargout=delete(dgObj,varargin)
%
% DELETE deletes specific results stored in OCOBJ.
%
% DELETE(OCOBJ,FIELDNAME,INDEX) deletes the results stored under FIELDNAME
% in the Result filed of OCOBJ. INDEX is an array specifying the index of
% the results to be deleted. 
%
% DELETE(OCOBJ,FIELDNAME1,INDEX1,FIELDNAME2,INDEX2, ...)
%
% OCOBJN=DELETE(...) the result is returned to the new stdocmodel instance
% OCOBJN

if isempty(dgObj)
    return
end
resultStruct=result(dgObj);
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
dgObj.Result=resultStruct;


if ~nargout
    assignin('caller',inputname(1),dgObj);
else
    varargout{1}=dgObj;
end

end
