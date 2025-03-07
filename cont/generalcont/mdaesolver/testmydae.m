function test=testmydae(flag,varargin)
%
% for the computation the solution components are reordered such that the
% first order components are at the beginnig and thereafter the zero order
% components. The sol.y components are in the initial ordering.
% If sol.xnew exist, then the solution has to be considered on the new
% meshgrid sol.xnew
test=[];
global OCMATCONT

switch flag
    case 'coeff2point'
        if nargin<2
            return
        end
        tcolmesh=varargin{1};
        coeff=varargin{2};
        [tcolmesh1_2,coeff1_2]=coeff2coeff1_2(tcolmesh,coeff);
        test.Ycf=coeff2points(tcolmesh,coeff,'collocationgrid2collocationgrid1_2');
        test.Ycf2=coeff2points(tcolmesh,coeff,'general',tcolmesh1_2);
        
        test.Yfc=coeff2points(tcolmesh1_2,coeff1_2,'collocationgrid1_22collocationgrid');
        test.Yfc2=coeff2points(tcolmesh1_2,coeff1_2,'general1_2',tcolmesh);

        order=[ones(OCMATCONT.firstordercomponentnumber,1);zeros(OCMATCONT.zeroordercomponentnumber,1)];
        order=order([OCMATCONT.firstordercoordinate;OCMATCONT.zeroordercoordinate]);
        tic
        test.YYcf=coeffToValues(coeff,tcolmesh(OCMATCONT.MESHDATA.tmeshidx),order,OCMATCONT.CollocationPoint,tcolmesh1_2);
        tst=toc;
        fprintf('Time duration for ''coarse to fine'': %f\n',tst)

        tic
        test.YYfc=coeffToValues(coeff1_2,tcolmesh1_2(OCMATCONT.MESHDATA1_2.tmeshidx),order,OCMATCONT.CollocationPoint,tcolmesh);
        tst=toc;
        fprintf('Time duration for ''fine to coarse'': %f\n',tst)
  
    case 'grid_grid1_2'
        tcolmesh=varargin{1};
        coeff=varargin{2};
        [tcolmesh1_2,coeff1_2]=coeff2coeff1_2(tcolmesh,coeff);
        test.y0=coeff2points(tcolmesh,coeff,'grid');
        test.y0col=coeff2points(tcolmesh,coeff,'collocationgrid');

        test.y=coeff2points(tcolmesh1_2,coeff1_2,'grid1_2');
        test.ycol=coeff2points(tcolmesh1_2,coeff1_2,'collocationgrid1_2');
        test.ycol2=coeff2points(tcolmesh1_2,coeff1_2,'general1_2',tcolmesh1_2);

        order=[ones(OCMATCONT.firstordercomponentnumber,1);zeros(OCMATCONT.zeroordercomponentnumber,1)];
        order=order([OCMATCONT.firstordercoordinate;OCMATCONT.zeroordercoordinate]);
        
        test.yy0=coeffToValues(coeff,tcolmesh(OCMATCONT.MESHDATA.tmeshidx),order,OCMATCONT.CollocationPoint,tcolmesh(OCMATCONT.MESHDATA.tmeshidx));
        test.yy0col=coeffToValues(coeff,tcolmesh(OCMATCONT.MESHDATA.tmeshidx),order,OCMATCONT.CollocationPoint,tcolmesh);
        
        test.yy=coeffToValues(coeff1_2,tcolmesh1_2(OCMATCONT.MESHDATA1_2.tmeshidx),order,OCMATCONT.CollocationPoint,tcolmesh1_2(OCMATCONT.MESHDATA1_2.tmeshidx));
        test.yycol=coeffToValues(coeff1_2,tcolmesh1_2(OCMATCONT.MESHDATA1_2.tmeshidx),order,OCMATCONT.CollocationPoint,tcolmesh1_2);
    case 'coeff2coeff1_2'
        order=[ones(OCMATCONT.firstordercomponentnumber,1);zeros(OCMATCONT.zeroordercomponentnumber,1)];
        order=order([OCMATCONT.firstordercoordinate;OCMATCONT.zeroordercoordinate]);

        tcolmesh=varargin{1};
        coeff=varargin{2};
        [tcolmesh1_2,coeff1_2]=coeff2coeff1_2(tcolmesh,coeff);
        test.coeff1_2=coeff1_2;
        
        
        initProfile.initialMesh=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);
        initProfile.initialValues=coeffToValues(coeff,tcolmesh(OCMATCONT.MESHDATA.tmeshidx),order,OCMATCONT.CollocationPoint);
        test.initProfile=initial_coefficients('bvps2_test_probdef',makefinemesh(tcolmesh(OCMATCONT.MESHDATA.tmeshidx)),initProfile,OCMATCONT.CollocationPoint,1);

    case 'point2coeff'
        if nargin<3
            return
        end
        sol=varargin{1};
        tnew=varargin{2};
        
    case 'collocationgrid'
        coeff=varargin{1};

        test.yfineT=coeff2points(coeff,'collocationgrid',[],1);
        test.yfine0=coeff2points(coeff,'collocationgrid_',[],1);

        coeff1_2=coeff2coeff1_2(coeff);
        coeff1_2=rand(length(coeff1_2),1);
        test.yfineT1_2=coeff2points(coeff1_2,'collocationgrid1_2',[],1);
        test.yfine01_2=coeff2points(coeff1_2,'collocationgrid1_2_',[],1);
        
    case 'collocationgridmT'
        coeff=varargin{1};
        test.yfine=coeff2points(coeff,'collocationgrid',[],1);
        test.yfinemT=coeff2points(coeff,'collocationgridmT',[],1);

    case 'midpoint'
        tcol=varargin{1};
        n=length(tcol);
        tic
        test.tcolmid = (tcol(1:end-1)+tcol(2:end))/2;
        tst=toc;fprintf('Time duration: %f\n',tst)
        tic
        test.tcolmid = tcol(1:end-1)+diff(tcol)/2;
        tst=toc;
        fprintf('Time duration: %f\n',tst)
        tic
        test.tcolmid = (tcol(1:n-1)+tcol(2:n))/2;
        tst=toc;
        fprintf('Time duration: %f\n',tst)
        tic
        test.tcolmid = tcol(1:n-1)+diff(tcol)/2;
        tst=toc;
        fprintf('Time duration: %f\n',tst)

    case 'residual'
        coeff=varargin{2};
        tcolmesh=varargin{1};
        funch=extremaldae4ft();
        tres=(tcolmesh(1:end-1)+tcolmesh(2:end))/2;
        test.res=funch{12}(tcolmesh,coeff,[],tres,funch{4}{1});
        test.resbvp=computeResidual('bvps2_mom_probdef',coeff,[],tcolmesh(OCMATCONT.MESHDATA.tmeshidx),OCMATCONT.CollocationPoint,tres,1);
        
                
        test.y=coeff2points(tcolmesh,coeff,'general',tres);

        order=[ones(OCMATCONT.firstordercomponentnumber,1);zeros(OCMATCONT.zeroordercomponentnumber,coeff(end))];
        order=order([OCMATCONT.firstordercoordinate;OCMATCONT.zeroordercoordinate]);
        test.yy=coeffToValues(coeff,tcolmesh(OCMATCONT.MESHDATA.tmeshidx),order,OCMATCONT.CollocationPoint,tres);

    case 'residual2'
        coeff=varargin{2};
        tmesh=varargin{1};
        tcolmesh=makecollocationmesh(tmesh);
        tres=(tcolmesh(1:end-1)+tcolmesh(2:end))/2;
             
        tic
        test.y_global=coeff2points_global(tcolmesh,coeff,'general',tres);
        tst=toc;
        fprintf('Time duration: %f\n',tst)
        tic
        test.y_new=coeff2points_new(tmesh,coeff,'general',tres);
        tst=toc;
        fprintf('Time duration: %f\n',tst)
             
        tic
        test.y2_global=coeff2points_global(tcolmesh,coeff,'grid');
        tst=toc;
        fprintf('Time duration: %f\n',tst)
        tic
        test.y2_new=coeff2points_new(tmesh,coeff,'grid');
        tst=toc;
        fprintf('Time duration: %f\n',tst)
             
        tic
        test.y3_global=coeff2points_global(tcolmesh,coeff,'collocationgrid');
        tst=toc;
        fprintf('Time duration: %f\n',tst)
        tic
        test.y3_new=coeff2points_new(tmesh,coeff,'collocationgrid');
        tst=toc;
        fprintf('Time duration: %f\n',tst)

end