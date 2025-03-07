function fn=filename(mmObj,varargin)
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
% the IDX'th parameter variables and values are taken for the file
% composition.
%
% FN = GETFILENAME(OCOBJ,FORMAT,IDX,PHIDX) PHIDX is an integer valued
% vector. Instead of the parameter value the placeholder operator * is used
% for the PHIDX'th variable.

fn='';
format='';
idx=[];
placeholderidx=[];

maxnumpar=3;
nummod=numberofmodels(mmObj);

if isempty(mmObj)
    return
end

if nargin>=2
   format=varargin{1};
end
if nargin>=3
   idx=varargin{2};
   if ischar(idx) || iscell(idx)
       idx=parameterindex(mmObj,idx);
   end
end
if nargin>=4
   placeholderidx=varargin{3};
end
if ischar(placeholderidx)
    placeholderidx=parameterindex(mmObj,placeholderidx);
end
if iscell(placeholderidx)
    placeholderidxtmp=placeholderidx;
    placeholderidx=[];
    for ii=1:numel(placeholderidxtmp)
        placeholderidx(ii)=parameterindex(mmObj,placeholderidxtmp{ii});
    end
end
if isempty(format)
    format='%3.4g';
end

% declare constants
[par,parvar]=parametervalue(mmObj);
fn=modelname(mmObj);
if isempty(idx)
    idx=repmat({1:maxnumpar},1,nummod);
end
if ischar(format)
    format=repmat(cellstr(format),1,length(idx));
end
if length(format)==1
    format=repmat(format,1,length(idx));
elseif length(format)~=length(idx)
    ocmaterror('Length of ''format'' cell differs from length of index.')
end

for kk=1:nummod
    fn=[fn '_P-' num2str(kk)];
    ii=0;
    for jj=idx{kk}
        ii=ii+1;
        if isempty(intersect(placeholderidx,jj))
            fn=[fn '_' parvar{kk}{jj} '_' num2str(par{kk}(jj),format{kk})];
        else
            fn=[fn '_' parvar{kk}{jj} '_*'];
        end
    end
end