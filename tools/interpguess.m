function [tmesh coeff]=interpguess(sol,bvpmethod)
%INTERP_GUESS  Evaluate/interpolate the initial guess at collocation points

global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

switch bvpmethod
    case 'bvp5c'
        sq5=sqrt(5);
        c=[ (5 - sq5)/10, (5 + sq5)/10];
        neqn=size(sol.y,1);

        h=diff(sol.x);
        x0 =sol.x(1:end-1);  % left subinerval boundaries
        xc1=x0 + c(1)*h;
        xc2=x0 + c(2)*h;
        try
            param=sol.parameters(:);
        catch
            param=[];
        end
        tmesh=[x0; xc1; xc2];
        tmesh=[tmesh(:); sol.x(end)];
        tmesh=reshape(tmesh,1,[]);

        y0=sol.y(:,1:end-1);  % solution at left subinterval boundaries
        if ~isempty(sol.solver)
            bvpmethod=sol.solver;
        else
            bvpmethod='unknown';
        end

        switch bvpmethod
            case 'unknown'
                % linear interpolation between mesh points
                dely=diff(sol.y,1,2);
                yc1=y0 + c(1)*dely;
                yc2=y0 + c(2)*dely;
            case 'bvp4c'
                % evaluate solution in SOL
                yc1=deval_bvp4c(sol.x,sol.y,xc1,sol.solverinfo.yp);
                yc2=deval_bvp4c(sol.x,sol.y,xc2,sol.solverinfo.yp);
            case 'bvp5c'
                % evaluate solution in SOL
                yc1=deval_bvp5c(sol.x,sol.y,xc1,sol.solverinfo.yp,sol.solverinfo.ymid);
                yc2=deval_bvp5c(sol.x,sol.y,xc2,sol.solverinfo.yp,sol.solverinfo.ymid);
            case 'bvp6c'
                % evaluate solution in SOL
                yc1=deval_bvp6c(sol.x,sol.y,xc1,sol.solverinfo.yp,sol.solverinfo.ypmid);
                yc2=deval_bvp6c(sol.x,sol.y,xc2,sol.solverinfo.yp,sol.solverinfo.ypmid);
        end
        yend=sol.y(:,end);

        coeff=[y0; yc1; yc2];
        coeff=[coeff(:); yend];
        coeff=[reshape(coeff,neqn,[]);param(:,ones(1,numel(tmesh)))];
        coeff=coeff(:);
    case 'sbvpoc'
        nummeshintv=OCMATCONT.HE.TIMEDDATA.nummeshintv;
        diffmesh=OCMATCONT.HE.TIMEDDATA.diffmesh;
        leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
        rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
        coeff=zeros(OCMATCONT.HE.numdvariables,1);
        counter=0;
        for ii=1:OCMATCONT.HE.numarc
            arcmesh=sol.x(leftarcindex(ii):rightarcindex(ii));
            arcdiffmesh=diffmesh(leftarcindex(ii):rightarcindex(ii)-1);
            arcindex=arcarg2arcindex(sol.arcarg(ii));
            numcols=domainddata(arcindex).numcols;
            eqcoord=domainddata(arcindex).eqcoord;
            odecoord=domainddata(arcindex).odecoord;
            aecoord=domainddata(arcindex).aecoord;
            numeq=domainddata(arcindex).numeq;
            numode=domainddata(arcindex).numode;
            numae=domainddata(arcindex).numae;
            numcolscoord=domainddata(arcindex).numcolscoord;
            psival=domainddata(arcindex).psivalequi;
            psival0=domainddata(arcindex).psivalequi0;

            if ~numae % no algebraic equation (implicit control) and ODEs are of order 1
                equidistmesh=[0:1/(numcols):1].';
                interpolationnummesh=arcmesh(ones(numcols,1),1:end-1)+equidistmesh(1:numcols,ones(1,nummeshintv(ii))).*arcdiffmesh(ones(numcols,1),:);
                interpolationnummesh=[interpolationnummesh(:).' arcmesh(end)];
                interpolatedval=interp1(arcmesh,sol.y(eqcoord,leftarcindex(ii):rightarcindex(ii)).',interpolationnummesh,'spline').';
                A=zeros(numcols+1);
                H=ones(numcols+1);
                A(1:numcols+1,1)=1; % coefficient for y_ijk
                A(:,numcolscoord+1)=reshape(psival(1,numcolscoord,1:numcols+1),numcols,numcols+1).'; % coefficient for z_ijlk
                for jj=0:nummeshintv(ii)-1
                    H(2:numcols+1,2:numcols+1)=arcdiffmesh(jj+1);
                    intervaldynVar=interpolatedval(:,jj*(numcols)+1:(jj+1)*(numcols)+1);
                    ((A.*H)\intervaldynVar').';
                    coeff([(jj*(numeq)*(numcols+1))+1:(jj+1)*((numeq)*(numcols+1))]+counter,1)=ans(:);
                end
                counter=(jj+1)*((numeq)*(numcols+1))+counter;


            else
                % for ODEs and AEs
                equidistmesh=[0:1/(numcols):1].';
                interpolationmesh=arcmesh(ones(numcols,1),1:end-1)+equidistmesh(1:numcols,ones(1,nummeshintv(ii))).*arcdiffmesh(ones(numcols,1),:);
                interpolationmesh=[interpolationmesh(:).' arcmesh(end)];
                interpolatedval=interp1(arcmesh,sol.y(odecoord,leftarcindex(ii):rightarcindex(ii)).',interpolationmesh,'spline').';
                equidistmeshae=[0:1/(numcols-1):1].';
                interpolationmeshae=arcmesh(ones(numcols-1,1),1:end-1)+equidistmeshae(1:numcols-1,ones(1,nummeshintv(ii))).*arcdiffmesh(ones(numcols-1,1),:);
                interpolationmeshae=[interpolationmeshae(:).' arcmesh(end)];
                if numae==1
                    interpolatedvalae=interp1(arcmesh,sol.y(aecoord,leftarcindex(ii):rightarcindex(ii)).',interpolationmeshae,'spline');
                else
                    interpolatedvalae=interp1(arcmesh,sol.y(aecoord,leftarcindex(ii):rightarcindex(ii)).',interpolationmeshae,'spline').';
                end
                Aae=reshape(psival0(1,numcolscoord,1:numcols),numcols,numcols).';
                A=zeros(numcols+1);
                H=ones(numcols+1);
                A(1:numcols+1,1)=1; % coefficient for y_ijk
                A(:,numcolscoord+1)=reshape(psival(1,numcolscoord,1:numcols+1),numcols,numcols+1).'; % coefficient for z_ijlk
                for jj=0:nummeshintv(ii)-1
                    counter_start=counter+1;
                    counter=counter+numode*(numcols+1)+numae*numcols;
                    H(2:numcols+1,2:numcols+1)=arcdiffmesh(jj+1);
                    intervaldynVar=interpolatedval(odecoord,jj*(numcols)+1:(jj+1)*(numcols)+1);
                    Code=((A.*H)\intervaldynVar').';
                    intervaldynVarae=interpolatedvalae(1:numae,jj*(numcols-1)+1:(jj+1)*(numcols-1)+1);
                    Cae=((Aae)\intervaldynVarae').';
                    C=[Code(odecoord,2:numcols+1);Cae];
                    coeff(counter_start:counter,1)=[Code(odecoord,1);C(:)];
                end
            end

        end
        coeff(OCMATCONT.HE.parametercoord,1)=sol.parameters(:);
        tmesh=sol.x;
end