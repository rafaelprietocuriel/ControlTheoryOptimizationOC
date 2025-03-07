function initocmatconthybridvariable(sol,codimension,contcoordinate)
global OCMATCONT 

if nargin==2
    codimension=1;
end
% discretization data of the domain
domainddata=OCMATCONT.DOMAINDDATA;
switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        % initialize global data for the time mesh
        OCMATCONT.HE.TIMEDDATA.mesh=sol.x;
        OCMATCONT.HE.TIMEDDATA.nummesh=numel(sol.x);
        OCMATCONT.HE.TIMEDDATA.diffmesh=diff(sol.x);

        % determine the number of arcs and the corresponding time mesh,
        % switches of arcs are indicated by a doubling of the switch time
        arcposition=find(OCMATCONT.HE.TIMEDDATA.diffmesh==0);
        OCMATCONT.HE.TIMEDDATA.leftarcindex=[1 arcposition+1];
        OCMATCONT.HE.TIMEDDATA.rightarcindex=[arcposition OCMATCONT.HE.TIMEDDATA.nummesh];
        OCMATCONT.HE.TIMEDDATA.nummeshintv=OCMATCONT.HE.TIMEDDATA.rightarcindex-OCMATCONT.HE.TIMEDDATA.leftarcindex;
        OCMATCONT.HE.numarc=numel(arcposition)+1;

        OCMATCONT.HE.arcindex=arcarg2arcindex(sol.arcarg);
        OCMATCONT.HE.arcarg=sol.arcarg;
        % edge [arcarg1 arcarg2] denotes the switch from arc with
        % arc argument arcarg1 to arc argument arcarg2, depending on this
        % information e.g. the boundary condtions have to be chosen
        OCMATCONT.HE.edge=[OCMATCONT.HE.arcarg(1:end-1).' OCMATCONT.HE.arcarg(2:end).'];
        OCMATCONT.HE.numparameter=numel(sol.parameters);%additional/continuation parameters
        OCMATCONT.HE.numparametermc=OCMATCONT.HE.numparameter-codimension;%additional parameters minus continuation parameter

        % determine the total number of unknowns and calculate the
        % collocation points
        numdvariables=0;

        % to make the following formulas more clear we redefine some of the
        % variables which in principle is superfluous
        nummeshintv=OCMATCONT.HE.TIMEDDATA.nummeshintv;
        tmesh=OCMATCONT.HE.TIMEDDATA.mesh;
        diffmesh=OCMATCONT.HE.TIMEDDATA.diffmesh;
        leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
        rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
        totalnumode=0;
        totalnumeq=0;
        maxnumeq=-inf;
        arccoeffindex=ones(OCMATCONT.HE.numarc,2);
        for arc=1:OCMATCONT.HE.numarc
            arcindex=OCMATCONT.HE.arcindex(arc);
            numcols=domainddata(arcindex).numcols;
            numeq=domainddata(arcindex).numeq;
            numode=domainddata(arcindex).numode;
            nodes=domainddata(arcindex).nodes;
            numcolscoord=domainddata(arcindex).numcolscoord;
            totalnumode=totalnumode+numode;
            totalnumeq=totalnumeq+numeq;
            maxnumeq=max(maxnumeq,numeq);
            % numeq*numcols at each collocation point the ODEs and AEs have
            % to be satisfied
            % at the mesh points (for each arc) the solution of the ODEs
            % are continuously connected
            numdvariables=numdvariables+(numeq*numcols+numode)*nummeshintv(arc);
            % the coordinates of the variables corresponding to the arc ARC
            % arccoeffindex(ARC,1):arccoeffindex(ARC,2)
            arccoeffindex(arc,2)=numdvariables;
            if arc<OCMATCONT.HE.numarc
                arccoeffindex(arc+1,1)=arccoeffindex(arc,2)+1;
            end
            OCMATCONT.HE.DDATA(arc).collocationpoints=tmesh(leftarcindex(arc):rightarcindex(arc)-1).'*ones(1,numcols)+diffmesh(leftarcindex(arc):rightarcindex(arc)-1).'*nodes(numcolscoord);
        end
        OCMATCONT.HE.arccoeffindex=arccoeffindex;
        OCMATCONT.HE.numdvariables=numdvariables+OCMATCONT.HE.numparameter;
        OCMATCONT.HE.totalnumode=totalnumode; % total number of equations for ODEs (does not include algebraic equations)
        OCMATCONT.HE.totalnumeq=totalnumeq; % total number of equations for ODEs (does not include algebraic equations)
        OCMATCONT.HE.maxnumeq=maxnumeq; % total number of equations for ODEs and AEs

        % HE.DDATA(ARC).meshvalcoord ... coordinates for the mesh points of
        % arc ARC
        % HE.DDATA(ARC).collvalcoord ... coordinates for the collocation
        % points of arc ARC
        % HE.parametercoord ... coordinates of the free and continuation
        % parameter, free parameter are e.g. unknown switching times
        % HE.contparametercoord ... coordinates of the continuation
        % parameter
        % HE.coeffcoordmp ... last coordinate excluding the free and
        % continuation parameter
        counter=0;
        IDX=1:OCMATCONT.HE.numdvariables;
        for arc=1:OCMATCONT.HE.numarc
            arcindex=OCMATCONT.HE.arcindex(arc);
            numcols=domainddata(arcindex).numcols;
            numeq=domainddata(arcindex).numeq;
            eqcoord=domainddata(arcindex).eqcoord;
            numcolscoord=domainddata(arcindex).numcolscoord;
            for jj=1:nummeshintv(arc)
                for ll=eqcoord
                    if domainddata(arcindex).order(ll)==1
                        counter=counter+1;
                        OCMATCONT.HE.DDATA(arc).meshvalcoord(ll,1,jj)=IDX(counter);
                    end
                end
                add2counter=domainddata(arcindex).numeq*domainddata(arcindex).numcols;
                OCMATCONT.HE.DDATA(arc).collvalcoord(eqcoord,numcolscoord,jj)=reshape(IDX(counter+1:counter+add2counter),numeq,numcols);
                counter=counter+add2counter;
            end
        end
        OCMATCONT.HE.parametercoord=IDX(counter+1:counter+OCMATCONT.HE.numparameter).';
        OCMATCONT.HE.contparametercoord=IDX(counter+OCMATCONT.HE.numparameter-1+(1:codimension)).';
        OCMATCONT.HE.coeffcoordmp=OCMATCONT.HE.numdvariables-OCMATCONT.HE.numparameter; % coordinates of coefficients without parameter values

        [b,msg]=testconsistency(OCMATCONT,sol);
        if ~b
            ocmaterror(msg)
        end

    case {'bvp4c','bvp5c','bvp6c'}

end
