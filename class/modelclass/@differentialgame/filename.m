function fn=filename(dgObj,varargin)
%
% GETFILENAME returns the filename of the oc model
%
% FN = GETFILENAME(OCOBJ) returns the filname composed of the name of the
% model and its parameters.
%
% FN = GETFILENAME(OCOBJ,FORMAT) uses the format string FORMAT (see SPRINTF
% for details) for the parameter values.
%
% FN = GETFILENAME(OCOBJ,FORMAT,IDX) IDX is an integer valued vector. Only
% the IDX'th parameter variables and values aretaken forthe file
% composition.
%
% FN = GETFILENAME(OCOBJ,FORMAT,IDX,PHIDX) PHIDX is an integer valued
% vector. Instead of the parameter value the placeholder operator * is used
% for the PHIDX'th variable.

fn='';
format='';
idx=[];
placeholderidx=[];
if isempty(dgObj)
    return
end

if nargin>=2
   format=varargin{1};
end
if nargin>=3
   idx=varargin{2};
   if ischar(idx)
       idx=parameterindex(dgObj,idx);
   end
end
if nargin>=4
   placeholderidx=varargin{3};
end
if ischar(placeholderidx)
    placeholderidx=parameterindex(dgObj,placeholderidx);
end
if iscell(placeholderidx)
    placeholderidxtmp=placeholderidx;
    placeholderidx=[];
    for ii=1:numel(placeholderidxtmp)
        placeholderidx(ii)=parameterindex(dgObj,placeholderidxtmp{ii});
    end
end
if isempty(format)
    format='%3.4g';
end

% declare constants
[par parvar]=parametervalue(dgObj);
fn=modelname(dgObj);
if isempty(idx)
    idx=1:numel(par);
end
if ischar(format)
    format=repmat(cellstr(format),1,length(idx));
end
if length(format)==1
    format=repmat(format,1,length(idx));
elseif length(format)~=length(idx)
    ocmaterror('Length of ''format'' cell differs from length of index.')
end

ii=0;
for jj=idx
    ii=ii+1;
    if isempty(intersect(placeholderidx,jj))
        fn=[fn '_' parvar{jj} '_' num2str(par(jj),format{ii})];
    else
        fn=[fn '_' parvar{jj} '_*'];
    end
end