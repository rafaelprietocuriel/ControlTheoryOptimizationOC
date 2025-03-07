function solnew=devalgrad(sol,tnew,xnew,method1d,method2d,varargin)
%

solnew=sol;
if nargin==3
    method1d='linear';
end
if isempty(method1d)
    method1d='linear';
end
if nargin<=4
    method2d=method1d;
end
if isempty(method2d)
    method2d=method1d;
end
tinterpflag=false;
if ischar(tnew)
    fac=str2double(tnew(1:end-1));
    tnew=linspace(sol.t(1),sol.t(end),fac*length(sol.t));
end
if ischar(xnew)
    fac=str2double(xnew(1:end-1));
    xnew=linspace(sol.x(1),sol.x(end),fac*length(sol.x));
end
if ~isempty(tnew)
    try
        if any(sol.t-tnew)
            tinterpflag=true;
        end
    catch
        tinterpflag=true;
    end
end

xinterpflag=false;
if ~isempty(xnew)
    try
        if any(sol.x-xnew)
            xinterpflag=true;
        end
    catch
        xinterpflag=true;
    end
end

if ~tinterpflag && ~xinterpflag
    solnew=sol;
    return
end

solnew=rmfield(solnew,'numericalinfo');
solnew=rmfield(solnew,'solverinfo');
if all([tinterpflag xinterpflag])
    [T,X]=meshgrid(sol.t,sol.x);
    solnew.t=tnew;
    solnew.x=xnew;
    nt=length(tnew);
    nx=length(xnew);
    solnew.u=zeros(nx,size(sol.u,2),nt);
    for ii=1:size(sol.u,2)
        solnew.u(:,ii,:)=reshape(interp2(T,X,squeeze(sol.u(:,ii,:)),solnew.t(:).',solnew.x(:),method2d),[nx,1,nt]);
    end
    solnew.y=zeros(nx,size(sol.y,2),nt);
    for ii=1:size(sol.y,2)
        solnew.y(:,ii,:)=reshape(interp2(T,X,squeeze(sol.y(:,ii,:)),solnew.t(:).',solnew.x(:),method2d),[nx,1,nt]);
    end
    if ~isempty(sol.v)
        solnew.v=interp1(sol.t,sol.v.',solnew.t,method1d).';
    end
    if ~isempty(sol.z)
        solnew.z=interp1(sol.t,sol.z.',solnew.t,method1d).';
    end
else
    if tinterpflag
        solnew.t=tnew;
        solnew.u=zeros(size(sol.u,1),size(sol.u,2),nt);
        for ii=1:size(sol.u,2)
            solnew.u(:,ii,:)=interp1(sol.t,squeeze(sol.u(:,ii,:)).',solnew.t,method1d).';
        end
        solnew.y=zeros(size(sol.y,1),size(sol.y,2),nt);
        for ii=1:size(sol.y,2)
            solnew.y(:,ii,:)=interp1(sol.t,squeeze(sol.y(:,ii,:)).',solnew.t,method1d).';
        end
        if ~isempty(sol.v)
            solnew.v=interp1(sol.t,sol.v.',solnew.t,method1d).';
        end
        if ~isempty(sol.z)
            solnew.z=interp1(sol.t,sol.z.',solnew.t,method1d).';
        end
    else
        solnew.x=xnew;
        nx=length(xnew);
        solnew.u=zeros(nx,size(sol.u,2),size(sol.u,3));
        for ii=1:size(sol.u,2)
            solnew.u(:,ii,:)=interp1(sol.x,squeeze(sol.u(:,ii,:)),solnew.x,method1d);
        end
        solnew.y=zeros(nx,size(sol.y,2),size(sol.y,3));
        for ii=1:size(sol.y,2)
            solnew.y(:,ii,:)=interp1(sol.x,squeeze(sol.y(:,ii,:)),solnew.x,method1d);
        end
    end
end
if ~isempty(solnew.y)
    solnew.modelinfo.y0=solnew.y(:,1:end/2,1);
end
if ~isempty(solnew.z)
    solnew.modelinfo.z0=solnew.z(1:end/2,1);
end
