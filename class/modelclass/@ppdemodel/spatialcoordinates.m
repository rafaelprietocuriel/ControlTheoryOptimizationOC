function out=spatialcoordinates(ppdeObj,solObj)

out=[];
if isempty(ppdeObj)
    return
end
% coordinates of the spatial state and costate
femdat=femdata(solObj);
statenum=retrieveppdemodelinformation(ppdeObj.Model,'statenum');
costatenum=retrieveppdemodelinformation(ppdeObj.Model,'costatenum');
gridnum=femdat.gridnum;
spatialstatedependenceindex=retrieveppdemodelinformation(ppdeObj.Model,'spatialstatedependenceindex');
nonspatialstatedependenceindex=retrieveppdemodelinformation(ppdeObj.Model,'nonspatialstatedependenceindex');
spatialcostatedependenceindex=retrieveppdemodelinformation(ppdeObj.Model,'spatialcostatedependenceindex');

ctr=0;
base=0;
for ii=1:statenum.value
    if any(spatialstatedependenceindex.value==ii)
        out(1:gridnum,ii)=(ii-1)*gridnum+1:ii*gridnum;
    elseif any(nonspatialstatedependenceindex.value==ii)
        if ctr==0
            base=(ii-1)*gridnum;
        end
        ctr=ctr+1;
        
        base=base+1;
    end
end
if ctr==0
    if ii==1
        base=gridnum;
    else
        base=ii*gridnum;
    end
end
offset=ii-ctr;
for ii=1:costatenum.value
    if any(spatialcostatedependenceindex.value==ii)
        out(1:gridnum,offset+ii)=base+((ii-1)*gridnum+1:ii*gridnum);
    end
end