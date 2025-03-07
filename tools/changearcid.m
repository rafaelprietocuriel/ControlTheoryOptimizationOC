function sol=changearcid(sol,changearcidx,newarcid,varargin)
%
%

if isempty(sol) || isempty(changearcidx)
    return
end
sol.arcarg(changearcidx)=newarcid;
