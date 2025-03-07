function fn=filename(odeObj,varargin)
%
% GETFILENAME returns the filename of the oc model
%
% FN = GETFILENAME(ODEOBJ) returns the filname composed of the name of the
% model and its parameters.
%
% FN = GETFILENAME(ODEOBJ,FORMAT) uses the format string FORMAT (see SPRINTF for
% details) for the parameter values.

fn='';
format='';
idx=[];
placeholderidx=[];
if isempty(odeObj)
    return
end

if nargin>=2
   format=varargin{1};
end
if nargin>=3
   idx=varargin{2};
   if ischar(idx)
       idx=parameterindex(odeObj,idx);
   end
end
if nargin>=4
   placeholderidx=varargin{3};
end
if ischar(placeholderidx)
    placeholderidx=parameterindex(odeObj,placeholderidx);
end
if iscell(placeholderidx)
    placeholderidxtmp=placeholderidx;
    placeholderidx=[];
    for ii=1:numel(placeholderidxtmp)
        placeholderidx(ii)=parameterindex(odeObj,placeholderidxtmp{ii});
    end
end
if isempty(format)
    format='%3.4g';
end

% declare constants
[par parvar]=parametervalue(odeObj);
fn=modelname(odeObj);
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