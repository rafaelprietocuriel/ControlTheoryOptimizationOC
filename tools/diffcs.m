function [A,z]=diffcs(fun,t,dynVar,pararg,arcarg,type)
% DIFFCS    Complex Step Jacobian
% J = jacobiancsd(f,x) returns the numerical (m x n) Jacobian matrix of a
% m-vector function, f(x) at the reference point, x (n-vector).
% [J,z] = jacobiancsd(f,x) also returns the function value, z=f(x).
%
% Example
% f=@(x)[x(1)*x(2)^2;x(1)*x(3)^2;x(1)^2];
% x=1:3;
% [J,z]=jacobiancsd(f,x)
% J =
%     4     4     0
%     9     0     6
%     2     0     0
% z =
%     4
%     9
%     1
%
% based on jacobiancsd
% By Yi Cao at Cranfield University, 02/01/2008
%
z=fun(t,dynVar,pararg,arcarg);                       % function value
%z=fun(t,dynVar,pararg);                       % function value

switch type
    case 'state'
        n=numel(dynVar);                     % size of independent
        m=numel(z);                     % size of dependent
        A=zeros(m,n);                   % allocate memory for the Jacobian matrix
        h=n*eps;                        % differentiation step size
        for k=1:n                       % loop for each independent variable
            dynVar1=dynVar;                       % reference point
            dynVar1(k)=dynVar1(k)+h*i;            % increment in kth independent variable
            A(:,k)=imag(fun(t,dynVar1,pararg,arcarg))/h;     % complex step differentiation
            %A(:,k)=imag(fun(t,dynVar1,pararg))/h;     % complex step differentiation
        end
    case 'parameter'
        n=numel(pararg);                     % size of independent
        m=numel(z);                     % size of dependent
        A=zeros(m,n);                   % allocate memory for the Jacobian matrix
        h=n*eps;                        % differentiation step size
        for k=1:n                       % loop for each independent variable
            pararg1=pararg;                       % reference point
            pararg1(k)=pararg1(k)+h*i;            % increment in kth independent variable
            A(:,k)=imag(fun(t,dynVar,pararg1,arcarg))/h;     % complex step differentiation
        end
    case 'time'
        m=numel(z);                     % size of dependent
        A=zeros(m,1);                   % allocate memory for the Jacobian matrix
        h=eps;                        % differentiation step size
        t1=t;                       % reference point
        t1=t1+h*i;            % increment in kth independent variable
        A=imag(fun(t1,dynVar,pararg,arcarg))/h;     % complex step differentiation
end