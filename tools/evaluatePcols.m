function ret=evaluatePcols(ii,diffmesh,y,z,coord,psival,numcols)
%Polynomial optimized version

% evaluate Runge-Kutta-basis polynomial at collocation points
% z denotes the coefficients for the basis polynomial (integrated Lagrange
% polynomials), y denotes the coefficients for the constant part (in
% general t^(k-1)/(k-1)!, k=1,...,max order) 

% numcols ... number of collocation points
if ~isempty(y)
    % differential equation
    ret=y(coord,1,ii(ones(1,numcols)))+sum(diffmesh((ii))*z(coord,:,ii(ones(1,numcols))).*psival(ones(numel(coord),1),:,1:numcols),2);
    ret=ret(:,:);
else
    % algebraic equation
    ret=sum(z(coord,:,ii(ones(1,numcols))).*psival(ones(numel(coord),1),:,1:numcols),2);
    ret=ret(:,:);
end