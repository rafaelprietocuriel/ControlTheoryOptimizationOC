function rs=reachableset(ocObj,solObj,varargin)
%
% REACHABLESET
rs=[];
if isempty(ocObj)
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if isscalar(solObj)
    contRes=contresult(ocObj,solObj);
    [contSol,contClass]=contsolution(contRes{1});
else
    return
end

if isempty(contSol)
    return
else
    for ii=1:length(contSol)
        switch contClass
            case {'extremal4ft','extremalt4ft','limitextremal4ft'}
                ocSol=octrajectory(contSol(ii));
                X=state(ocObj,ocSol,1);
                rsStruct.y(:,ii)=X(:,end);
                rsStruct.userinfo.timehorizon(ii)=timehorizon(ocSol);
                if ii==1
                    rsStruct.arcarg=[];
                    rsStruct.arcposition=[];
                    rsStruct.userinfo.solutiontype='octrajectory';
                    rsStruct.userinfo.parametervalue=par;
                end
            case {'indifferencesolution4ft'}
                ocSol=sol2ocmultipath(contSol(ii));
                for jj=1:ocSol(1).solverinfo.indifferenceorder
                    X=state(ocObj,ocSol(jj),1);
                    rsStruct(jj).y(:,ii)=X(:,end);
                    rsStruct(jj).userinfo.timehorizon(ii)=timehorizon(ocSol(jj));
                    if ii==1
                        rsStruct(jj).arcarg=[];
                        rsStruct(jj).arcposition=[];
                        rsStruct(jj).userinfo.solutiontype='octrajectory';
                        rsStruct(jj).userinfo.parametervalue=par;
                    end
                end
        end
    end
end

for ii=1:length(rsStruct)
    rs{ii}=occurve(rsStruct(ii));
end

