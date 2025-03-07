function val=evalcollocationpoly(idx,coord,t,tmesh,y,z,diffmesh,numcolscoord,psi,order) %Polynomial of derivation abl

if order
    val=y(coord,idx);
    for ii=numcolscoord
        val=val+diffmesh(idx)*z(coord,ii,idx).*repmat(polyval(psi(ii,:),(t-tmesh(idx))/diffmesh(idx)),numel(coord),1);
    end
else
    val=0;
    for ii=numcolscoord
        val=val+z(coord,ii,idx).*repmat(polyval(psi(ii,:),(t-tmesh(idx))/diffmesh(idx)),numel(coord),1);
    end
end
