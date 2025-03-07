function femdata=generatefemdata(ppdeObj,gridnum,varargin)

femdata=[];
if isempty(ppdeObj)
    return
end

spacetype=spacegeometrytype(ppdeObj);

if strcmp(spacetype,'interval')
    oopdeobj=eval([modelname(ppdeObj) 'OOPDE1D']);
    par=parametervalue(ppdeObj);
    femdata.gridnum=gridnum;
    femdata.gridnumm1=gridnum-1;
    initialize(oopdeobj,par,femdata);
    femdata.K=oopdeobj.A;
    femdata.M=oopdeobj.M;
    femdata.invMK=femdata.M\femdata.K;
    femdata.grid=oopdeobj.grid;
    initialize(oopdeobj,par,(femdata.grid.x(1:end-1)+femdata.grid.x(2:end))*0.5);
    femdata.gridmid=oopdeobj.grid;
end