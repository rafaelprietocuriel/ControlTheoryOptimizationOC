function out=ellipticpde()
%
                                                                                                  
out{1}=@epde;
out{2}=@boundarycondition;
out{3}=@jacobian;
out{4}=@jacobianp;
out{5}=[];
out{6}=[];
out{7}=[];
out{8}=[];
out{9}=[];
out{10}= [];
out{20}= @plotcont;
%
% --------------------------------------------------------------------------
function res=epde(coeff,par,arcid)
global PPDEPRIMITIVE
[x,y]=rearr(coeff);
F=PPDEPRIMITIVE.pdefun([],x,y,par,arcid);
res=PPDEPRIMITIVE.femop.M*F-(PPDEPRIMITIVE.femop.K-PPDEPRIMITIVE.femop.Kadv)*coeff;

function out=boundarycondition(coeff,par,arcid)
global PPDEPRIMITIVE
out=PPDEPRIMITIVE.boundarycondition;

% --------------------------------------------------------------------------
function jac=jacobian(coeff,par,arcid)
global PPDEPRIMITIVE
[x,y]=rearr(coeff);

if PPDEPRIMITIVE.symjac % analytically
    Fu=PPDEPRIMITIVE.pdejacobian([],x,y,par,arcid);
    jac=PPDEPRIMITIVE.femop.M*Fu-(PPDEPRIMITIVE.femop.K-PPDEPRIMITIVE.femop.Kadv);
end

% --------------------------------------------------------------------------
function jacp=jacobianp(coeff,varargin)
global PPDEPRIMITIVE
[x,y]=rearr(coeff);
jacp=-PPDEPRIMITIVE.pdeparameterjacobian([],x,y,par,arcid);

function h=plotcont(coeff)
global PPDEPRIMITIVE
[x,y]=rearr(coeff);
h=[];
%h=plot(x(cds.ndim,1:contnum),x(PPDEPRIMITIVE.plotcoord,1:contnum));

function [x y]=rearr(coeff)
global PPDEPRIMITIVE
y=coeff(PPDEPRIMITIVE.coeffidx);
x=PPDEPRIMITIVE.points;