function out=trianglearea(ppdePrim)

pt=points(ppdePrim);
tr=triangles(ppdePrim);
out=pdetrg(pt,tr);
