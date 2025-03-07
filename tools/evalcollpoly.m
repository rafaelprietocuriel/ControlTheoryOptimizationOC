function ret=evalcollpoly(ii,diffmesh,y,z,arc,posflag,aeflag)
% simplified for order one autonomous ODEs

% evaluation of the "Runge-Kutta basis polynomial" at time one (end time of
% the ith interval). For ODEs this value is given by y. For algebraic
% equations it is the value of the Legendre polynomial.
global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;
arcindex=OCMATCONT.HE.arcindex(arc);
aecoord=domainddata(arcindex).aecoord;
odecoord=domainddata(arcindex).odecoord;
numae=domainddata(arcindex).numae;
numode=domainddata(arcindex).numode;
numcols=domainddata(arcindex).numcols;

switch posflag
    case 0
        if ~aeflag
            ret=y(odecoord,1,ii);
        else
            ret=sum(z(aecoord,:,ii).*domainddata(arcindex).psival0(ones(numae,1),:,numcols+2),2);
        end
    case 1
        if ~aeflag
            ret=y(odecoord,1,ii)+sum(diffmesh(ii+arc-1)*z(odecoord,:,ii).*domainddata(arcindex).psival(ones(numode,1),:,numcols+1),2);
        else
            ret=sum(z(aecoord,:,ii).*domainddata(arcindex).psival0(ones(numae,1),:,numcols+1),2);
        end
    case 0.5
end
