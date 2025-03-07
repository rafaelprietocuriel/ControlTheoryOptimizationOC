function ocgTrj=distance(ocgTrj1,ocgTrj2,varargin)

spec=[];
dim=[];
if nargin>=3
    spec=varargin{1};
end
if nargin>=4
    dim=varargin{2};
end
if isempty(spec)
    spec=1;
end
arcn1=arcnum(ocgTrj1);
arcn2=arcnum(ocgTrj2);
if arcn1~=arcn2
    ocgTrj=ocgtrajectory();
    return
end
if isempty(dim)
    dim=odenumber(ocgTrj1);
end
switch spec
    case 0
        ocgTrj=ocgTrj1;
        ocgTrj.octrajectory.solverinfo=[];
        for ii=1:arcn1
            ocgTrj.octrajectory.y{ii}=ocgTrj1.octrajectory.y{ii}(1:dim(ii),:)-ocgTrj2.octrajectory.y{ii}(1:dim(ii),:);
        end
    case inf
    otherwise
        ocgTrj=ocgTrj1;
        ocgTrj.octrajectory.solverinfo=[];
        for ii=1:arcn1
            ocgTrj.octrajectory.y{ii}=sum(abs(ocgTrj1.octrajectory.y{ii}(1:dim(ii),:)-ocgTrj2.octrajectory.y{ii}(1:dim(ii),:)).^spec).^(1/spec);
        end
end
