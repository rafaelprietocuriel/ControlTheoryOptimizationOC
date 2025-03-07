function out=cartel2dsym4IndifferenceSolutionContinuation()
%
% main function for the ocmat continuation process
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
out{1}=@canonicalsystem;
out{2}{1}=@jacobian;
out{2}{2}=@parameterjacobian;
	
out{3}{1}=@hessian;
out{3}{2}=@parameterhessian;
	
out{5}{1}=@bcinitial;
out{5}{2}=@bcasymptotic;
out{5}{3}=@bctransversality;
out{5}{4}=@equilibrium;
out{5}{5}=@bcindifference;
out{5}{6}=@salvagevalue;
out{5}{7}=@algebraicequation;
	
out{6}{1}=@bcjacobianinitial;
out{6}{2}=@bcjacobianasymptotic;
out{6}{3}=@bcjacobiantransversality;
out{6}{4}=@bcjacobianindifference;
	
% functions for hybrid structure and its discretization
out{7}{1}=@hybridinfo;
out{7}{2}=@domain;
out{7}{3}=@guard;
out{7}{4}=@reset;
out{7}{5}=@switchtime;
out{7}{7}=@jacobianguard;
out{7}{8}=@jacobianreset;
out{7}{9}=@domaindiscretization;
out{7}{10}=@timesettransformation;
out{8}{1}=@objectivefunction;
out{8}{2}=@objectivefunctionjacobian;
out{8}{3}=@objectivefunctionparameterjacobian;
out{8}{4}=@objectivefunctionderivativetime;

out{10}{1}=@targetfunction;
	
out{11}=@plotcontinuation;
out{12}=@testadmissibility;
	
out{20}=@datapath;
out{21}=@saveintermediatefiles;
	
function dxdt=canonicalsystem(t,depvar,par,arcid)
dxdt=cartel2dsymCanonicalSystemGeneralized(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function J=jacobian(t,depvar,par,arcid)
J=cartel2dsymCanonicalSystemGeneralizedJacobian(t,depvar,par,arcid);
	
	
%--------------------------------------------------------------------------
function J=parameterjacobian(t,depvar,par,arcid)
J=cartel2dsymCanonicalSystemParameterJacobian(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function H=hessian(t,depvar,par,arcid)
H=cartel2dsymCanonicalSystemHessian(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function H=parameterhessian(t,depvar,par,arcid)
H=cartel2dsymCanonicalSystemTotalHessian(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function O=objectivefunction(t,depvar,par,arcid)
O=cartel2dsymObjectiveFunction(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function J=objectivefunctionjacobian(t,depvar,par,arcid)
J=cartel2dsymObjectiveFunctionJacobian(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function Jpar=objectivefunctionparameterjacobian(t,depvar,par,arcid)
Jpar=cartel2dsymObjectiveFunctionParameterJacobian(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function Jt=objectivefunctionderivativetime(t,depvar,par,arcid)
Jt=cartel2dsymObjectiveFunctionDerivativeTime(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function res=equilibrium(depvar,par,arcid)
res=cartel2dsymEquilibriumEquation(depvar,par,arcid);
	
%--------------------------------------------------------------------------
function val=guard(depvara,depvarb,par,switchtimes,arcid,edge,arc)
% condition on the (continuous) variables (states, costates) when a
% switching from arcid1 to arcid2 takes place. edge=[arcid1 arcid2]
% this condition determines the switching time
if isempty(edge)
	val=[];
else
	val=cartel2dsymConstraint(switchtimes(arc),depvara(:,arc+1),par,edge(2,arc))-cartel2dsymConstraint(switchtimes(arc),depvarb(:,arc),par,edge(1,arc));
	% the user has to adapt the following lines to the constraint properties of the model
	if edge(1,arc)==0 && edge(2,arc)==1 || edge(2,arc)==0  && edge(1,arc)==1 || ...
            edge(1,arc)==3 && edge(2,arc)==6 || edge(2,arc)==3  && edge(1,arc)==6 || ...
            edge(1,arc)==4 && edge(2,arc)==7 || edge(2,arc)==4  && edge(1,arc)==7 || ...
            edge(1,arc)==8 && edge(2,arc)==11 || edge(2,arc)==8  && edge(1,arc)==11
		val=val(1);
	elseif edge(1,arc)==0 && edge(2,arc)==2 || edge(2,arc)==0  && edge(1,arc)==2 || ...
            edge(1,arc)==3 && edge(2,arc)==8 || edge(2,arc)==3  && edge(1,arc)==8  || ...
            edge(1,arc)==3 && edge(2,arc)==11 || edge(2,arc)==3  && edge(1,arc)==11    || ...
            edge(1,arc)==6 && edge(2,arc)==11 || edge(2,arc)==6  && edge(1,arc)==11     
		val=val(2);
    elseif  edge(1,arc)==4 && edge(2,arc)==10 || edge(2,arc)==4  && edge(1,arc)==10 || ...
            edge(1,arc)==0 && edge(2,arc)==3 || edge(2,arc)==0  && edge(1,arc)==3 || ...
            edge(1,arc)==1 && edge(2,arc)==6 || edge(2,arc)==1  && edge(1,arc)==6 || ...
            edge(1,arc)==7 && edge(2,arc)==13 || edge(2,arc)==7  && edge(1,arc)==13
        val=val(3);
	elseif edge(1,arc)==0 && edge(2,arc)==4 || edge(2,arc)==0  && edge(1,arc)==4 || ...
            edge(1,arc)==3 && edge(2,arc)==10 || edge(2,arc)==3  && edge(1,arc)==10 || ...
            edge(1,arc)==1 && edge(2,arc)==7 || edge(2,arc)==1  && edge(1,arc)==7 || ...
            edge(1,arc)==6 && edge(2,arc)==13 || edge(2,arc)==6  && edge(1,arc)==13
		val=val(4);
	end
end
	
%--------------------------------------------------------------------------
function val=reset(depvara,depvarb,par,switchtimes,arcid,edge,arc)
%transformation of the (continuous) variables (states, costates) when
%switching from arcid1 to arcid2. edge=[arcid1 arcid2] 
if isempty(edge)
	val=[];
else
	val=[depvara(1:4,arc+1)-depvarb(1:4,arc); ...
		cartel2dsymAlgebraicEquation(0,depvara(:,arc+1),par,edge(2,arc))];
end
%--------------------------------------------------------------------------
function val=targetfunction(depvara,depvarb,par,arcid)
% condition on the (continuous) variables (states, costates) at a switch to different stages
val=depvara(1,1)-depvara(2,1);
%--------------------------------------------------------------------------
function val=bcconnectingparts(depvara,depvarb,par,switchtimes,arcid,arc)
% condition on the (continuous) variables (states, costates) at a switch to different stages
val=depvara(1,1)-depvara(2,1);
%--------------------------------------------------------------------------
function val=bcoptimalconnectingparts(depvara,depvarb,par,arctim,arcid,arc)
% condition on the (continuous) variables (states, costates) at a switch to different stages that is optimal
val=cartel2dsymHamiltonian(arctim(arc+1),depvara(:,arc+1),par{arc+1},arcid{arc+1}(1))-cartel2dsymHamiltonian(arctim(arc+1),depvarb(:,arc),par{arc},arcid{arc}(end));
	
%--------------------------------------------------------------------------
function res=bcinitial(depvara,targetcoordinate,initialstate,par,arcid)
res=[depvara(targetcoordinate,1)-initialstate; ...
	cartel2dsymAlgebraicEquation(0,depvara(:,1),par,arcid)];
	
%--------------------------------------------------------------------------
function res=algebraicequation(depvara,par,arcid)
res=cartel2dsymAlgebraicEquation(0,depvara(:,1),par,arcid);
	
%--------------------------------------------------------------------------
function res=bcasymptotic(depvarb,asymptoticmatrix,saddlepoint)
res=asymptoticmatrix'*(depvarb(1:4,end)-saddlepoint(1:4));
	
%--------------------------------------------------------------------------
function res=bctransversality(T,depvarb,par,arcid)
res=depvarb(3:4,end)-cartel2dsymTransversalityBC(T,depvarb(:,end),par,arcid);
	
%--------------------------------------------------------------------------
function out=salvagevalue(T,depvarb,par,arcid)
out=cartel2dsymDiscountedSalvagevalue(T,depvarb,par,arcid);
	
%--------------------------------------------------------------------------
function res=bcindifference(depvara,par,arcid,initcoordinate)
res=cartel2dsymHamiltonian(0,depvara(:,initcoordinate(1)),par,arcid(initcoordinate(1)))-cartel2dsymHamiltonian(0,depvara(:,initcoordinate(2)),par,arcid(initcoordinate(2)));
	
%--------------------------------------------------------------------------
function out=objectivevalue(t,depvar,par,arcid)
out=cartel2dsymObjectivevalue(t,depvar,par,arcid);
	
%--------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjacobianinitial(depvara,freepar,arcid,targetcoordinate,continuationvector)
[Ja,Jb,Jpar]=cartel2dsymBCJacobian4Initial(depvara,freepar,arcid,targetcoordinate,continuationvector);
	
%--------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjacobianasymptotic(depvarb,arcid,asymptoticmatrix,saddlepoint)
[Ja,Jb,Jpar]=cartel2dsymBCJacobian4Asymptotic(depvarb,arcid,asymptoticmatrix,saddlepoint);
	
%--------------------------------------------------------------------------
function [Ja,Jb,Jpar]=bcjacobianindifference(depvara,par,arcid,initcoordinate)
[Ja,Jb,Jpar]=cartel2dsymBCJacobian4Indifference(depvara,par,arcid,initcoordinate);
	
%--------------------------------------------------------------------------
function J=jacobianguard(depvara,depvarb,par,switchtimes,arcid,edge,arc,arcoffset)
J=cartel2dsymBCJacobian4Guard(depvara,depvarb,par,switchtimes(arc),arcid,edge,arc,arcoffset);
	
%--------------------------------------------------------------------------
function J=jacobianreset(depvara,depvarb,par,switchtimes,arcid,edge,arc,arcoffset)
J=cartel2dsymBCJacobian4Reset(depvara,depvarb,par,switchtimes(arc),arcid,edge,arc,arcoffset);
	
%--------------------------------------------------------------------------
function infoH=hybridinfo()
infoH.arcarg=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]; % replaces name of arcnum arcidentifier
infoH.edge=[0 1;0 2;0 3;0 4;0 5;0 6;0 7;0 8;0 9;0 10;0 11;0 12;0 13;0 14;1 2;1 3;1 4;1 5;1 6;1 7;1 8;1 9;1 10;1 11;1 12;1 13;1 14;2 3;2 4;2 5;2 6;2 7;2 8;2 9;2 10;2 11;2 12;2 13;2 14;3 4;3 5;3 6;3 7;3 8;3 9;3 10;3 11;3 12;3 13;3 14;4 5;4 6;4 7;4 8;4 9;4 10;4 11;4 12;4 13;4 14;5 6;5 7;5 8;5 9;5 10;5 11;5 12;5 13;5 14;6 7;6 8;6 9;6 10;6 11;6 12;6 13;6 14;7 8;7 9;7 10;7 11;7 12;7 13;7 14;8 9;8 10;8 11;8 12;8 13;8 14;9 10;9 11;9 12;9 13;9 14;10 11;10 12;10 13;10 14;11 12;11 13;11 14;12 13;12 14;13 14;1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0;13 0;14 0;2 1;3 1;4 1;5 1;6 1;7 1;8 1;9 1;10 1;11 1;12 1;13 1;14 1;3 2;4 2;5 2;6 2;7 2;8 2;9 2;10 2;11 2;12 2;13 2;14 2;4 3;5 3;6 3;7 3;8 3;9 3;10 3;11 3;12 3;13 3;14 3;5 4;6 4;7 4;8 4;9 4;10 4;11 4;12 4;13 4;14 4;6 5;7 5;8 5;9 5;10 5;11 5;12 5;13 5;14 5;7 6;8 6;9 6;10 6;11 6;12 6;13 6;14 6;8 7;9 7;10 7;11 7;12 7;13 7;14 7;9 8;10 8;11 8;12 8;13 8;14 8;10 9;11 9;12 9;13 9;14 9;11 10;12 10;13 10;14 10;12 11;13 11;14 11;13 12;14 12;14 13]; % each row represents a possible transition between the states of the hybrid system    
	
%--------------------------------------------------------------------------
function infoD=domain(arcid)
% returns general information about the domain for each arc
infoD=cartel2dsymArcInfo(arcid);
	
%--------------------------------------------------------------------------
function val=switchtime(edge)
% here a fixed switchtime for a specific edge can be specified, otherwise
% it is determined by the guard
val=[];
	
%--------------------------------------------------------------------------
function infoD=domaindiscretization(arcid)
% returns information about the discretization for each arc
infoD=cartel2dsymArcDiscretizationInfo(arcid);
	
%--------------------------------------------------------------------------
function infoT=timesettransformation()
infoT.normalization=1; % the timenterval of the nth arc is normalized to [n,n+1]
%infoT.infinity=-0.1; % for an infinite time horizon problem the last interval [N,inf)->[N,N+1]
%infoT.asymptoticapproximation=inf;
infoT.infinity=0; % for an infinite time horizon problem the last interval [N,inf)->[N,N+1], infoT.infinity<0
infoT.asymptoticapproximation=100;
	
%--------------------------------------------------------------------------
function h=plotcontinuation(t,depvar,modelpar,arcid,freepar,tangent)
h=cartel2dsymPlotIndifferenceContinuation(t,depvar,modelpar,arcid,freepar,tangent);
	
%--------------------------------------------------------------------------
function [out labelS]=testadmissibility(t,depvar,par,arcid)
	
out=cartel2dsymAdmissible(t,depvar,par,arcid);
	
labelS(1).arcid=0;
labelS(1).info='';
labelS(2).arcid=1;
labelS(2).info='';
labelS(3).arcid=2;
labelS(3).info='';
labelS(4).arcid=3;
labelS(4).info='';
labelS(5).arcid=4;
labelS(5).info='';
labelS(6).arcid=5;
labelS(6).info='';
labelS(7).arcid=6;
labelS(7).info='';
labelS(8).arcid=7;
labelS(8).info='';
labelS(9).arcid=8;
labelS(9).info='';
labelS(10).arcid=9;
labelS(10).info='';
labelS(11).arcid=10;
labelS(11).info='';
labelS(12).arcid=11;
labelS(12).info='';
labelS(13).arcid=12;
labelS(13).info='';
labelS(14).arcid=13;
labelS(14).info='';
labelS(15).arcid=14;
labelS(15).info='';
	
%--------------------------------------------------------------------------
function pathname=datapath()
	
pathname=getocmatfolder('userdata','stdocmodel','cartel2dsym');
	
%--------------------------------------------------------------------------
function [resultfile,globalvarfile]=saveintermediatefiles()
	
resultfile='SaveIntermediateResults';
globalvarfile='SaveIntermediateResultsGlobalVariable';
