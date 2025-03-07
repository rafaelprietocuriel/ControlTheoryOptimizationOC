function idx=contsol2idx(contSol,fieldname,paridx,parval)
%
% CONT2IDX returns indices where a curve crosses a specific parameter value
%
% IDX=CONT2IDX(CONTPAR,PARVAL) returns the index/indices where the
% curve stored as a vector of discrete numbers CONTPAR, usually coming from
% a continuation process, crosses the parameter value PARVAL.

idx=[];

if isempty(parval)
            return
end

if ~isfield(contSol(1),fieldname) && ~isfield(contSol(1).solverinfo,fieldname)
    return
end

if isempty(paridx)
    paridx=1;
end
contpar=zeros(1,length(contSol));
for ii=1:length(contSol)
    if ~isfield(contSol(ii),fieldname)
        contpar(ii)=contSol(ii).solverinfo.(fieldname)(paridx);
    else
        contpar(ii)=contSol(ii).(fieldname)(paridx);
    end
end
if any(abs(imag(contpar))>0)
    warning('Parameter values have imaginary parts. ')
    contpar=real(contpar);
end
idx=find(diff(sign(contpar-parval)));
