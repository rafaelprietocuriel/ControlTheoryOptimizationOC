function tangentint=devaltangent(tmeshint,tmesh,tangent,neqn,N,neqnN)
global OCMATCONT OCBVP

if isempty(tangent)
    tangentint=[];
    return
end
mbcidx = find(diff(tmesh) == 0);  % locate internal interfaces
ismbvp = ~isempty(mbcidx);
Lidx = [1, mbcidx+1];
Ridx = [mbcidx, length(tmesh)];
nregions = length(mbcidx) + 1;
mbcidxint = find(diff(tmeshint) == 0);  % locate internal interfaces
Lidxint = [1, mbcidxint+1];
Ridxint = [mbcidxint, length(tmeshint)];
if length(neqn)>1
        tangentfull=zeros(OCBVP.maxnumode,length(tmesh));
        tangentfull(OCMATCONT.HE.DDATA.meshvalcoord)=tangent(OCMATCONT.HE.ycoord);
        tangentint=[];
        for region=1:nregions
            xidx=Lidx(region):Ridx(region);
            xidxint=Lidxint(region):Ridxint(region);
            tmp=interp1(tmesh(xidx),tangentfull(1:OCBVP.numode(region),xidx).',tmeshint(xidxint)).';
            tangentint=[tangentint;tmp(:)];
        end
        freepar=tangent(OCMATCONT.HE.parametercoord); % parameter values
else
    freepar=tangent(neqnN+1:end);
    tangent=reshape(tangent(1:neqnN),neqn,N);
    tangentint=[];
    if ismbvp
        for region=1:nregions
            xidx=Lidx(region):Ridx(region);
            xidxint=Lidxint(region):Ridxint(region);
            tangentint=[tangentint interp1(tmesh(xidx),tangent(:,xidx).',tmeshint(xidxint)).'];
        end
    else
        tangentint=interp1(tmesh,tangent.',tmeshint).';
    end
end
tangentint=[tangentint(:);freepar];
tangentint=tangentint/norm(tangentint);