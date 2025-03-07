function varargout=discretizationdata(tmesh,coeff,flag,varargin)
%
% CALC_DISCRETIZATIONDATA returns the data of the discretization
%
% flag is a string: 'x1tau' for the timemesh together with the collocation points
%           'valx1' returns the state values at the time mesh
%           'valx1tau' returns the state values at the time mesh and the
%           collocation points

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

switch OCMATCONT.solver
    case 'sbvpoc'
        [tmesh,y,z,par]=OCMATCONT.drearr(tmesh,coeff);
        diffmesh=diff(tmesh);

        switch flag
            case 'collmeshdata'
                % returns solution evaluated at mesh and collocation points
                % it has to be assured that the data sored in TIMEDDATA and
                % DDATA are consistent
                leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
                rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
                maxnumcols=max([domainddata(:).numcols]);
                removecolumn=[];
                collmesh=[];
                absshift=0;
                for arc=1:OCMATCONT.HE.numarc
                    arcindex=OCMATCONT.HE.arcindex(arc);
                    cols=OCMATCONT.HE.DDATA(arc).collocationpoints;
                    idx=leftarcindex(arc)-arc+1:rightarcindex(arc)-arc;
                    psival=domainddata(arcindex).psival;
                    psival0=domainddata(arcindex).psival0;
                    numcols=domainddata(arcindex).numcols;
                    diffnumcols=maxnumcols-numcols;
                    numcolscoord=domainddata(arcindex).numcolscoord;
                    numae=domainddata(arcindex).numae;
                    numode=domainddata(arcindex).numode;
                    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
                    aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations

                    if diffnumcols
                        counter=(maxnumcols+1).*(1:rightarcindex(arc)-leftarcindex(arc));
                        repmat(absshift,diffnumcols,rightarcindex(arc)-leftarcindex(arc))+repmat(counter,diffnumcols,1)+repmat((1:diffnumcols)'-diffnumcols,1,rightarcindex(arc)-leftarcindex(arc));
                        removecolumn=[removecolumn ans(:).'];
                    end
                    collmeshtmp=[tmesh(leftarcindex(arc):rightarcindex(arc)-1).' cols].';
                    collmesh=[collmesh collmeshtmp(:).' tmesh(rightarcindex(arc))];
                    collmeshval(odecoord,1,leftarcindex(arc):rightarcindex(arc)-1)=y(odecoord,1,idx);
                    if ~isempty(aecoord)
                        pos=numcols+2;
                        collmeshval(aecoord,1,leftarcindex(arc):rightarcindex(arc)-1)=sum(z(aecoord,numcolscoord,idx).*psival0(ones(numae,1),numcolscoord,pos(ones(1,rightarcindex(arc)-leftarcindex(arc)))),2);
                    end
                    for ii=numcolscoord
                        collmeshval(odecoord,ii+1,leftarcindex(arc):rightarcindex(arc)-1)=y(odecoord,1,idx)+reshape(diffmesh(ones(numode,1),leftarcindex(arc):rightarcindex(arc)-1),numode,1,rightarcindex(arc)-leftarcindex(arc)).*sum(z(odecoord,numcolscoord,idx).*psival(ones(numode,1),numcolscoord,ii(ones(1,rightarcindex(arc)-leftarcindex(arc)))),2);
                        if ~isempty(aecoord)
                            collmeshval(aecoord,ii+1,leftarcindex(arc):rightarcindex(arc)-1)=sum(z(aecoord,numcolscoord,idx).*psival0(ones(numae,1),numcolscoord,ii(ones(1,rightarcindex(arc)-leftarcindex(arc)))),2);
                        end
                    end
                    pos=numcols+1;
                    collmeshval(odecoord,1,rightarcindex(arc))=y(odecoord,1,idx(end))+diffmesh(rightarcindex(arc)-1)*sum(z(odecoord,numcolscoord,idx(end)).*psival(ones(numode,1),numcolscoord,pos),2);
                    if ~isempty(aecoord)
                        collmeshval(aecoord,1,rightarcindex(arc))=sum(z(aecoord,numcolscoord,idx(end)).*psival0(ones(numae,1),numcolscoord,pos),2);
                    end
                    removecolumn=[removecolumn (rightarcindex(arc)-1)*(maxnumcols+1)+1+[1:maxnumcols]];
                    absshift=removecolumn(end);

                end
                collmeshval=collmeshval(:,:);
                % for testing purposes
%                 if any(find(~sum(collmeshval))-removecolumn)
%                     warning('Syntax for removing zero xolumns is not correct')
%                 end
                collmeshval(:,removecolumn)=[];

                varargout{1}=collmesh;
                varargout{2}=collmeshval;
                varargout{3}=par;
            case 'colldata'
                % returns solution evaluated at mesh and collocation points
                % it has to be assured that the data sored in TIMEDDATA and
                % DDATA are consistent
                leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
                rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
                maxnumcols=max([domainddata(:).numcols]);
                removecolumn=[];
                absshift=0;
                for arc=1:OCMATCONT.HE.numarc
                    arcindex=OCMATCONT.HE.arcindex(arc);
                    cols=OCMATCONT.HE.DDATA(arc).collocationpoints;
                    idx=leftarcindex(arc)-arc+1:rightarcindex(arc)-arc;
                    psival=domainddata(arcindex).psival;
                    psival0=domainddata(arcindex).psival0;
                    numcols=domainddata(arcindex).numcols;
                    diffnumcols=maxnumcols-numcols;
                    numcolscoord=domainddata(arcindex).numcolscoord;
                    numae=domainddata(arcindex).numae;
                    numode=domainddata(arcindex).numode;
                    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
                    aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations

                    if diffnumcols
                        counter=(maxnumcols+1).*(1:rightarcindex(arc)-leftarcindex(arc));
                        repmat(absshift,diffnumcols,rightarcindex(arc)-leftarcindex(arc))+repmat(counter,diffnumcols,1)+repmat((1:diffnumcols)'-diffnumcols,1,rightarcindex(arc)-leftarcindex(arc));
                        removecolumn=[removecolumn ans(:).'];
                    end
                    collmesh=[tmesh(leftarcindex(arc):rightarcindex(arc)-1).' cols].';
                    collmesh=[collmesh(:).' tmesh(rightarcindex(arc))];
                    for ii=numcolscoord
                        collval(odecoord,ii,idx)=y(odecoord,1,idx)+reshape(diffmesh(ones(numode,1),leftarcindex(arc):rightarcindex(arc)-1),numode,1,rightarcindex(arc)-leftarcindex(arc)).*sum(z(odecoord,numcolscoord,idx).*psival(ones(numode,1),numcolscoord,ii(ones(1,rightarcindex(arc)-leftarcindex(arc)))),2);
                        if ~isempty(aecoord)
                            collval(aecoord,ii,idx)=sum(z(aecoord,numcolscoord,idx).*psival0(ones(numae,1),numcolscoord,ii(ones(1,rightarcindex(arc)-leftarcindex(arc)))),2);
                        end
                    end
                    absshift=(rightarcindex(arc)-arc)*maxnumcols;
                end
                collval=collval(:,:);
                % for testing purposes
%                 if any(find(~sum(collval))-removecolumn)
%                     warning('Syntax for removing zero columns is not correct')
%                 end
                collval(:,removecolumn)=[];
                varargout{1}=collmesh;
                varargout{2}=collval;
            case 'meshdata'
                leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
                rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
                for arc=1:OCMATCONT.HE.numarc
                    arcindex=OCMATCONT.HE.arcindex(arc);
                    idx=leftarcindex(arc)-arc+1:rightarcindex(arc)-arc;
                    psival=domainddata(arcindex).psival;
                    psival0=domainddata(arcindex).psival0;
                    numcols=domainddata(arcindex).numcols;
                    numcolscoord=domainddata(arcindex).numcolscoord;
                    numode=domainddata(arcindex).numode;
                    numae=domainddata(arcindex).numae;
                    odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
                    aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations

                    meshval(odecoord,1,leftarcindex(arc):rightarcindex(arc)-1)=y(odecoord,1,idx);
                    if ~isempty(aecoord)
                        pos=numcols+2;
                        meshval(aecoord,1,leftarcindex(arc):rightarcindex(arc)-1)=sum(z(aecoord,numcolscoord,idx).*psival0(ones(numae,1),numcolscoord,pos(ones(1,rightarcindex(arc)-leftarcindex(arc)))),2);
                    end
                    pos=numcols+1;
                    meshval(odecoord,1,rightarcindex(arc))=y(odecoord,1,idx(end))+diffmesh(rightarcindex(arc)-1)*sum(z(odecoord,numcolscoord,idx(end)).*psival(ones(numode,1),numcolscoord,pos),2);
                    if ~isempty(aecoord)
                        meshval(aecoord,1,rightarcindex(arc))=sum(z(aecoord,numcolscoord,idx(end)).*psival0(ones(numae,1),numcolscoord,pos),2);
                    end
                end
                varargout{1}=meshval(:,:);
        end
end
